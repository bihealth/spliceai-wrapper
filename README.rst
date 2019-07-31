================
SpliceAI Wrapper
================

.. image:: https://img.shields.io/pypi/v/spliceai-wrapper.svg
        :target: https://pypi.python.org/pypi/spliceai-wrapper

.. image:: https://img.shields.io/travis/bihealth/spliceai-wrapper.svg
        :target: https://travis-ci.org/bihealth/spliceai-wrapper

`Illumina SpliceAI <https://github.com/Illumina/SpliceAI>`_ is a nice method for predicting the impact of variants on splicing.
However, it is computationally very expensive (45k variants/hour on a GPU, a few hundred variants per hour and CPU core).

This project, **SpliceAI Wrapper**, is an attempt to use caching for reducing the number of required predictions.

------------
Installation
------------

I recommend to use Bioconda

.. code-block:: bash

    $ conda install spliceai-wrapper

If you're not installing from Bioconda, make sure that you have ``bcftools`` and ``spliceai`` installed and the executables in your path.

----------------------------
Importing Precomputed Scores
----------------------------

First, obtain the precomputed scores from the SpliceAI project (I'm using the genome-wide ones filtered to a score >= 0.1 for space usage reasons).
Then:

.. code-block:: bash

    $ spliceai-wrapper prepare \
        --release GRCh37 \
        --precomputed-db-path path/to/precomputed.sqlite3 \
        --precomputed-vcf-path path/to/whole_genome_filtered_spliceai_scores.vcf.gz

This will import the precomputed scores into a SQLite3 database.
On my workstation, it takes about 20 minutes.

------------------------
Running SpliceAI Wrapper
------------------------

Obtain the gene list text file from the SpliceAI project.
Then:

.. code-block:: bash

    $ spliceai-wrapper annotate \
        --input-vcf INPUT.vcf.gz \
        --output-vcf OUTPUT.vcf.gz \
        --genes-tsv path/to/grch37.txt \
        --precomputed-db-path path/to/precomputed.sqlite3 \
        --cache-db-path path/to/cache.sqlite3 \
        --path-reference path/to/hs37d5.fa \
        --release GRCh37

For trying it out use the ``--head 500`` parameter.

This will first go through ``INPUT.vcf.gz`` and try to find precomputed or cached values for all variants.
These precomputed/cached values will be used for annotation.
Variants that lie outside the genes defined in ``grch37.txt`` are ignored.
SNVs that lie within the genes defined in ``grch37.txt`` and that are not precomputed will be ignored as well (it is assumed their score is <0.1 otherwise they would appear).

The remaining variants will be written to a temporary VCF file and ``spliceai`` will be called on them.
The annotations from the output of ``spliceai`` will be cached and the output VCF file and VCF file with cache hits will be merged into ``OUTPUT.vcf.gz``.

Notes:

- The precomputation database is opened read-only so you also don't need write permissions to this file.
- The cache file must be writeable by your user, of course.
- The extension of your output file determines what format is used for writing it.
  ``.bcf`` files are written as compressed BCF, ``.vcf.gz`` and ``.vcf.bgz`` are written as bgzip-ed VCF, all other files will be written as text VCF files.
