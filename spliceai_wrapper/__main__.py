# -*- coding: utf-8 -*-
"""Console entry point for SpliceAI Wrapper."""

import argparse
from contextlib import closing, ExitStack
import os.path
import pathlib
import subprocess
import sqlite3
import sys
import tempfile
import vcfpy

from logzero import logger
from ncls import NCLS
import numpy as np
from tqdm import tqdm
from xdg import XDG_CACHE_HOME, XDG_DATA_HOME

from spliceai_wrapper import __version__


#: SQLite3 query for creating precomputation database.
#:
#: Indices are the same as in VCF, 1-based.
SQL_CREATE_TABLE = r"""
CREATE TABLE IF NOT EXISTS %(release)s_%(table_name)s
(
    var_desc TEXT PRIMARY KEY,

    chromosome VARCHAR(64),
    position INTEGER,
    reference TEXT,
    alternative TEXT,

    symbol TEXT,
    strand CHARACTER,
    var_type CHARACTER,
    distance INTEGER,
    delta_score_acceptor_gain FLOAT,
    delta_score_acceptor_loss FLOAT,
    delta_score_donor_gain FLOAT,
    delta_score_donor_loss FLOAT,
    delta_position_acceptor_gain INTEGER,
    delta_position_acceptor_loss INTEGER,
    delta_position_donor_gain INTEGER,
    delta_position_donor_loss INTEGER
);
"""

#: SQLite3 query for inserting data.
SQL_UPSERT = r"""
INSERT INTO %(release)s_%(table_name)s
(
    var_desc,
    chromosome,
    position,
    reference,
    alternative,
    symbol,
    strand,
    var_type,
    distance,
    delta_score_acceptor_gain,
    delta_score_acceptor_loss,
    delta_score_donor_gain,
    delta_score_donor_loss,
    delta_position_acceptor_gain,
    delta_position_acceptor_loss,
    delta_position_donor_gain,
    delta_position_donor_loss
)
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
ON CONFLICT DO NOTHING;
"""

#: SQLite3 query for selecting data.
SQL_SELECT = r"""
SELECT
    var_desc,
    alternative,
    symbol,
    delta_score_acceptor_gain,
    delta_score_acceptor_loss,
    delta_score_donor_gain,
    delta_score_donor_loss,
    delta_position_acceptor_gain,
    delta_position_acceptor_loss,
    delta_position_donor_gain,
    delta_position_donor_loss
FROM %(release)s_%(table_name)s
WHERE (var_desc IN (%%(var_descs)s));
"""


def yield_from_precomputed_vcf(reader):
    """Yield precomputed tuples from precomputed VCF file in ``reader``."""
    for record in reader:
        for alt in record.ALT:
            yield (
                "-".join(map(str, (record.CHROM, record.POS, record.REF, record.ALT[0].value))),
                record.CHROM,
                record.POS,
                record.REF,
                alt.value,
                record.INFO["SYMBOL"],
                record.INFO["STRAND"],
                record.INFO["TYPE"],
                record.INFO["DIST"],
                record.INFO["DS_AG"],
                record.INFO["DS_AL"],
                record.INFO["DS_DG"],
                record.INFO["DS_DL"],
                record.INFO["DP_AG"],
                record.INFO["DP_AL"],
                record.INFO["DP_DG"],
                record.INFO["DP_DL"],
            )


def yield_from_spliceai(reader):
    """Yield precomputed tuples from SpliceAI VCF in ``reader``."""
    for record in reader:
        for anno in record.INFO["SpliceAI"]:
            (alt, symbol, ds_ag, ds_al, ds_dg, ds_dl, dp_ag, dp_al, dp_dg, dp_dl) = anno.split("|")
            yield (
                "-".join(map(str, (record.CHROM, record.POS, record.REF, alt))),
                record.CHROM,
                record.POS,
                record.REF,
                alt,
                symbol,
                ".",  # don't know about strand
                ".",  # don't know about exonic/intronic
                -1,  # don't know about distance to splice site
                ds_ag,
                ds_al,
                ds_dg,
                ds_dl,
                dp_ag,
                dp_al,
                dp_dg,
                dp_dl,
            )


def prepare(args):
    logger.info("Running 'prepare' with args = %s", vars(args))

    if not os.path.exists(args.precomputed_db_path):
        os.makedirs(os.path.dirname(args.precomputed_db_path), exist_ok=True)

    logger.info("Opening database file %s", args.precomputed_db_path)
    with closing(sqlite3.connect(args.precomputed_db_path), isolation_level=None) as con:
        sql_kwargs = {"release": args.release, "table_name": "spliceai_scores"}
        sql_script = SQL_CREATE_TABLE % sql_kwargs
        logger.info("Executing %s to create table..." % sql_script)
        con.executescript(sql_script)

        logger.info("Opening VCF for import: %s...", args.precomputed_vcf_path)
        with vcfpy.Reader.from_path(args.precomputed_vcf_path) as reader:
            it = tqdm(yield_from_precomputed_vcf(reader), unit="records")
            con.executemany(SQL_UPSERT % sql_kwargs, it)
    logger.info("Done running 'prepare'.")


def augment_record(sql_records, vcf_record):
    records = {sql_record[0]: sql_record[1:] for sql_record in sql_records}
    keys = [
        "-".join(map(str, (vcf_record.CHROM, vcf_record.POS, vcf_record.REF, alt.value)))
        for alt in vcf_record.ALT
    ]
    arr = ["|".join(map(str, records[key])) for key in keys if key in records]
    if arr:
        vcf_record.INFO["SpliceAI"] = arr
    return vcf_record


def split_vcf(args, ncls, reader, con_pre, con_cache, writer_hit, writer_nohit):
    """Split VCF."""
    sql_select = SQL_SELECT % {"release": args.release, "table_name": "spliceai_scores"}
    cache_try = 0
    cache_hit = 0
    pre_try = 0
    pre_hit = 0
    pre_low = 0
    nohit = 0
    no_gene = 0
    for i, record in tqdm(enumerate(reader), unit="records"):
        if args.head and i >= args.head:
            break
        keys = [
            "-".join(map(str, (record.CHROM, record.POS, record.REF, alt.value)))
            for alt in record.ALT
        ]
        overlap_genes = record.CHROM in ncls and bool(
            list(ncls[record.CHROM].find_overlap(record.affected_start, record.affected_end))
        )
        is_snv = all(len(alt.value) == 1 for alt in record.ALT)
        # Query precomputed data.
        qry = sql_select % {"var_descs": ",".join(["?"] * len(keys))}
        if overlap_genes:
            if is_snv:  # precomputed only exists for in-gene SNVS
                pre_try += 1
                res_pre = list(con_pre.execute(qry, keys).fetchall())
                if res_pre:  # found precomputed
                    pre_hit += 1
                    writer_hit.write_record(augment_record(res_pre, record))
                else:  # not found, but must be cached => score <0.1
                    pre_low += 1
                    writer_hit.write_record(record)
                continue
            else:  # not SNV, but overlaps with genes, try to find in cache
                # Query cached data.
                cache_try += 1
                res_cache = list(con_cache.execute(qry, keys).fetchall())
                if res_cache:
                    cache_hit += 1
                    writer_hit.write_record(augment_record(res_cache, record))
                    continue
                # Write to no-hit file
                nohit += 1
                writer_nohit.write_record(record)
        else:
            # does not overlap with genes => write out to hit & ignore
            no_gene += 1
            writer_hit.write_record(record)
    logger.info(
        "Hits: %d/%d (%.1f%%), pre hits %d/%d (%.1f%%), pre low-score %d/%d (%.1f%%), cache "
        "hits %d/%d (%.1f%%), no gene: %d, cache misses: %d",
        cache_hit + pre_hit,
        cache_hit + pre_hit + pre_low + nohit,
        100.0 * (cache_hit + pre_hit) / (cache_hit + pre_hit + pre_low + nohit)
        if cache_hit + pre_hit + pre_low + nohit
        else 0.0,
        pre_low,
        cache_hit + pre_hit + pre_low + nohit,
        100.0 * pre_low / (cache_hit + pre_hit + pre_low + nohit)
        if cache_hit + pre_hit + pre_low + nohit
        else 0.0,
        pre_hit,
        pre_try,
        100.0 * pre_hit / pre_try if pre_try else 0.0,
        cache_hit,
        cache_try,
        100.0 * cache_hit / cache_try if cache_try else 0.0,
        no_gene,
        cache_try - cache_hit,
    )
    return cache_try - cache_hit  # == cache_misses


def annotate(args):
    logger.info("Loading gene intervals...")
    with open(args.genes_tsv) as inputf:
        intervals = {}
        header = None
        for i, line in enumerate(inputf):
            arr = line.strip().split("\t")
            if not header:
                arr[0] = arr[0][1:]
                header = arr
            else:
                record = dict(zip(header, arr))
                if record["CHROM"] not in intervals:
                    intervals[record["CHROM"]] = {"starts": [], "ends": [], "ids": []}
                entry = intervals[record["CHROM"]]
                entry["starts"].append(int(record["TX_START"]) - 1)
                entry["ends"].append(int(record["TX_END"]))
                entry["ids"].append(i)
    logger.info("Building NCLS...")
    ncls = {
        chrom: NCLS(
            np.asarray(entry["starts"]), np.asarray(entry["ends"]), np.asarray(entry["ids"])
        )
        for chrom, entry in intervals.items()
    }

    logger.info("Running 'annotate' with args = %s", vars(args))
    with ExitStack() as stack:
        logger.info("Opening %s (read-only)", args.precomputed_db_path)
        logger.info("URL = %s", "file:%s?mode=ro" % args.precomputed_db_path)
        con_pre = sqlite3.connect(
            "file:%s?mode=ro" % args.precomputed_db_path, uri=True, isolation_level=None
        )
        stack.enter_context(closing(con_pre))

        if not os.path.exists(args.cache_db_path):
            os.makedirs(os.path.dirname(args.cache_db_path), exist_ok=True)

        logger.info("Opening %s (cache; writeable)", args.cache_db_path)
        con_cache = sqlite3.connect(str(args.cache_db_path), isolation_level=None)
        stack.enter_context(closing(con_cache))
        sql_kwargs = {"release": args.release, "table_name": "spliceai_scores"}
        sql_script = SQL_CREATE_TABLE % sql_kwargs
        logger.info("Executing %s..." % sql_script)
        con_cache.executescript(sql_script)

        logger.info("Creating temporary directory...")
        tmpdir = tempfile.TemporaryDirectory()
        tmpdirname = stack.enter_context(tmpdir)
        logger.info(" => %s", tmpdirname)

        path_cache_hit = pathlib.Path(tmpdirname, "cache_hit.vcf")
        path_cache_nohit = pathlib.Path(tmpdirname, "cache_nohit.vcf")
        logger.info("Splitting %s", args.input_vcf)
        logger.info("  into: %s", path_cache_hit)
        logger.info("  and:  %s", path_cache_nohit)

        reader = stack.enter_context(vcfpy.Reader.from_path(args.input_vcf))
        header_hit = reader.header.copy()
        header_hit.add_info_line(
            {
                "ID": "SpliceAI",
                "Number": ".",
                "Type": "String",
                "Description": (
                    "SpliceAIv1.2.1 variant annotation. These include delta scores (DS) and delta positions (DP) "
                    "for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). "
                    "Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL"
                ),
            }
        )
        header_nohit = reader.header.copy()

        # Open vcfpy writers in a separate context, so everything is written and flushed once
        # split_vcf() is done.
        with ExitStack() as writer_stack:
            writer_hit = vcfpy.Writer.from_path(path_cache_hit, header_hit)
            writer_stack.enter_context(writer_hit)
            writer_nohit = vcfpy.Writer.from_path(path_cache_nohit, header_nohit)
            writer_stack.enter_context(writer_nohit)
            cache_misses = split_vcf(
                args, ncls, reader, con_pre, con_cache, writer_hit, writer_nohit
            )

        # Prepare arguments for bcftools
        out_fmt = "v"
        if args.output_vcf.endswith(".vcf.gz") or args.output_vcf.endswith(".vcf.bgz"):
            out_fmt = "z"
        elif args.output_vcf.endswith(".bcf"):
            out_fmt = "b"

        # Actually run spliceai.
        if not cache_misses:
            logger.info("No cache misses, no need to run spliceai")
            cmd = [
                "bcftools",
                "view",
                "-O",
                out_fmt,
                "-o",
                args.output_vcf,
                # the actual input files
                path_cache_hit,
            ]
            cmd = list(map(str, cmd))
            logger.info("Converting result with %s" % " ".join(cmd))
            p = subprocess.run(cmd)
            if p.returncode != 0:
                raise RuntimeError(
                    "Process exited with retcode %d: %s" % (p.returncode, " ".join(cmd))
                )
        else:
            path_spliceai_out = pathlib.Path(tmpdirname, "spliceai_out.vcf")
            cmd = [
                "spliceai",
                "-I",
                path_cache_nohit,
                "-O",
                path_spliceai_out,
                "-A",
                args.release.lower(),
                "-R",
                args.path_reference,
            ]
            cmd = list(map(str, cmd))
            logger.info("Running SpliceAI with %s" % " ".join(cmd))
            p = subprocess.run(cmd)
            if p.returncode != 0:
                raise RuntimeError(
                    "Process exited with retcode %d: %s" % (p.returncode, " ".join(cmd))
                )

            # Store results in cache.
            logger.info("Opening VCF for import: %s...", path_spliceai_out)
            with vcfpy.Reader.from_path(path_spliceai_out) as reader:
                it = tqdm(yield_from_spliceai(reader), unit="records")
                con_cache.executemany(SQL_UPSERT % sql_kwargs, it)

            # Use bcftools concat for merging the file to the output file.
            cmd = [
                "bcftools",
                "concat",
                "-O",
                out_fmt,
                "-o",
                args.output_vcf,
                # the actual input files
                path_cache_hit,
                path_spliceai_out,
            ]
            cmd = list(map(str, cmd))
            logger.info("Merging result with %s" % " ".join(cmd))
            p = subprocess.run(cmd)
            if p.returncode != 0:
                raise RuntimeError(
                    "Process exited with retcode %d: %s" % (p.returncode, " ".join(cmd))
                )

    logger.info("Done running 'annotate'.")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Caching wrapper for Illumina SpliceAI")
    parser.add_argument("--version", action="version", version="%%(prog)s %s" % __version__)
    subparsers = parser.add_subparsers(dest="action")

    parser_prepare = subparsers.add_parser(
        "prepare", help="Construct SQLite database from precomputed data"
    )
    parser_prepare.add_argument("--release", default="GRCh37", help="Release to use.")
    parser_prepare.add_argument(
        "--precomputed-db-path",
        default=(pathlib.Path(XDG_DATA_HOME) / "spliceai-wrapper" / "precomputed.sqlite3"),
    )
    parser_prepare.add_argument(
        "--precomputed-vcf-path",
        required=True,
        help="Path to VCF file for loading precomputed data from",
    )

    parser_annotate = subparsers.add_parser(
        "annotate", help="Annotate VCF file with SpliceAI using cache for the scores"
    )
    parser_annotate.add_argument(
        "--genes-tsv", required=True, help="Path to grch3[78].txt from SpliceAI"
    )
    parser_annotate.add_argument("--release", default="GRCh37", help="Release to use.")
    parser_annotate.add_argument(
        "--precomputed-db-path",
        default=(pathlib.Path(XDG_DATA_HOME) / "spliceai-wrapper" / "precomputed.sqlite3"),
    )
    parser_annotate.add_argument(
        "--cache-db-path",
        default=(pathlib.Path(XDG_CACHE_HOME) / "spliceai-wrapper" / "cache.sqlite3"),
        help="Path to SQLite3 file for the cache (to be updated)",
    )
    parser_annotate.add_argument("--input-vcf", required=True, help="Path to VCF file to annotate")
    parser_annotate.add_argument(
        "--output-vcf", required=True, help="Path to write annotated VCF to"
    )
    parser_annotate.add_argument(
        "--min-score",
        default=0.1,
        type=float,
        help="Minimal score to consider (report as 0 if smaller).",
    )
    parser_annotate.add_argument(
        "--head", default=None, type=int, help="Optional; only consider top N records."
    )
    parser_annotate.add_argument(
        "--path-reference", required=True, help="Path to reference FASTA file."
    )

    args = parser.parse_args(argv)

    if args.action == "prepare":
        return prepare(args)
    elif args.action == "annotate":
        return annotate(args)
    else:
        parser.print_help(sys.stderr)
        parser.exit("Invalid action: %s" % args.action)


if __name__ == "__main__":
    sys.exit(main())
