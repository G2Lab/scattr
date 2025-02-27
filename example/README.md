# Example of using ScatTR

## Running the example

In this example, we are estimating the copy number of a single TR in a simulated WGS sample. To run the example, make sure you have `scattr` or `cargo` installed in your `PATH`:

```
./run.sh
```

## Input files

In the `input/` directory, there are example input files to run ScatTR. In general, running ScatTR requires the following input files:

- A tab-delimeted **catalog** containing TR loci of interest: `catalog.tsv`
    - The file is expected to have columns: `id`, `contig`, `start`, `end`, and `motif`
- Aligned reads of the **sample** and their index: `input/sample.bam` and `input/sample.bam.bai`
- A **reference genome** and its index: `input/reference.fa` and `input/reference.fa.fai`

> **Note**: the simulated sample has only one contig (or chromosome) in its reference, namely `test_genome`. When running `scattr stats`, the option `--include-contigs test_genome` has to be added to correctly get the expected insert and depth distributions.

## Output files

In the `output/` directory, there are output files generated by ScatTR based on the example input files. ScatTR produces the following output files:

- The extracted bag of reads: `sample.bag.bam`
- The extracted depth and insert distributions: `sample.stats.json`
    - By default, ScatTR produces plots of these distributions (`sample.depth_distr.png` and `sample.insert_distr.png`). The `--no-plot` option disables this behavior
- The optimization problem definitions: `sample.defs.json`
    - This file information about all the reads associated with each TR and their relative positions with the decoy references (refer to manuscript methods section for more details)
- The estimated copy numbers: `sample.genotypes.json`
    - This file contains estimates for each TR locus in addition to a 95% confidence interval

## Generating input files

Generating the input files requires having `python`, `art_illumina`, `bwa`, `samtools` in your `PATH`. The input files are generated by running:

```
./scripts/make_inputs.sh
```

The script creates a random reference genome (~200 kbp in length) that is used to generate reads for a heteroyzygous TR expansion with normal allele copy number being 3 and with expanded copy number being 100. The motif of the TR is `CAGATA`.
