.. _getting-started:

===============
Getting started
===============

Setting up the environment
==========================

We use conda to ensure to create an isolated environment containing all of the
Python packages and external tools needed for the analyses. To create the conda
environment used for our analyses, run the following command from the sb-screen
root directory:

.. code-block:: bash

    make env

This creates an environment call *sb-screen*, which can be activated as follows:

.. code-block:: bash

    source activate sb-screen

Generating the datasets
=======================

Pre-computed
------------

To avoid re-generating the entire dataset (e.g. the insertion analysis, the
RNA-seq alignment and quantification), we provide a pre-computed version
of the dataset. This dataset is the same version as that was used for our
publication and should therefore give results identical or highly similar
to our analyses. The dataset can be downloaded as follows:

.. code-block:: bash

    rm -rf data/processed
    make download_freeze

Note that this command first clears the current data directory, before
downloading and extracting the frozen dataset.

From scratch
------------

To help generate the dataset from scratch, we provide a number of commands
in the Makefile that execute various scripts and Snakemake pipelines to
generate the different data types. Note that some steps of these pipelines
are stochastic (such as the CIS and NMF analyses). Therefore, results
obtained from the regenerated files may not be identical to our published
results. However, they should be highly similar.

Insertion analysis
~~~~~~~~~~~~~~~~~~

To run the insertion analysis, we first download the Ensembl mm10 reference
genome and use this reference to create an index for Bowtie2. In a second
step, we use PyIM (which requires the Bowtie2 index) to identify insertions
in the 99 samples that show an ILC morphology (see our publication for
more details). In a third step, we re-run this analysis using the full
set of 123 samples, which is required for an additional analysis. These three
steps are run using the following make commands:

.. code-block:: bash

    # Build the reference index
    make build_bowtie_index

    # Run the ShearSplink analysis in PyIM
    make shear_splink

    # Re-run the analysis using all samples
    make shear_splink_full

For the ShearSplink commands, extra options can be passed to the underlying
Snakemake pipeline using the SNAKEMAKE_ARGS variable. For example, to specify
the number of cores used by Snakemake when running the pipeline, you can use
the following command:

.. code-block:: bash

    make shear_splink SNAKEMAKE_ARGS='-j 5'

Similarly, a dry-run of the pipeline can be performed using:

.. code-block:: bash

    make shear_splink SNAKEMAKE_ARGS='-n -p'

RNA-seq analysis
~~~~~~~~~~~~~~~~

To run the RNA-seq analysis, we first download the reference genome (if this
was not already done before) and use the reference to generate an index for
STAR.

.. code-block:: bash

    make build_star_index

Next, to ensure that we can run the STAR efficiently in parallel, we
load this reference genome into shared memory before running the alignment.
After this, we can run the RNA-seq pipeline for the SB tumors using the
corresponding command. Together, this looks as follows:

.. code-block:: bash

    make preload_star_index
    make rnaseq_sb

Similarly, we can run the RNA-seq analysis for the KB1P and Pten mouse models
using the following commands:

.. code-block:: bash

    # Run the analysis for the KB1P samples
    make preload_star_index
    make rnaseq_kb1p

    # Run the analysis for the Pten samples
    make preload_star_index
    make rnaseq_pten

Note that we need to preload the index for each run of the RNA-seq pipeline,
as this pipeline currently unloads the shared reference during its execution.
This may be modified in a newer version of the pipeline. However, the extra
preloads should be instantaneous, as the shared reference is already in memory.

Similar to the insertion pipeline, the RNA-seq pipeline is a Snakemake pipeline
that can be passed extra arguments to enable parallel execution, etc.:

.. code-block:: bash

    make rnaseq_sb SNAKEMAKE_ARGS='-j 5'


NMF analysis
~~~~~~~~~~~~

Because of its computationally intensive nature, the NMF factorization used
in the subtype analysis is also pre-computed using an R script. This script
can be run after generating the SB RNA-seq dataset, using the following command:

.. code-block:: bash

    make nmf

The script performs both the rank estimation and the final factorization for
four clusters.
