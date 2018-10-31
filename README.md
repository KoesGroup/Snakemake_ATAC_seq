# ChIP_seq_Snakemake
A snakemake pipeline for the analysis of ATAC-seq data

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.2.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Miniconda](https://img.shields.io/badge/miniconda-blue.svg)](https://conda.io/miniconda)

# Aim

Snakemake pipeline made for reproducible analysis of paired-end Illumina ATAC-seq data. This pipeline is similar on many point to the Snakemake_ChIP-seq pipeline, it includes steps of mapping to the mitochondrial and chloroplast genome of tomato as contamination by these organelles occurs often in ATAC-seq.
Post-processing also includes the tool **AlignementSieve**, for the adjustement of the binding site.


For now I chose to include the fasta files for the mitochondrial and chloroplast genome directly in the repository, in the future I would like the pipeline to be able to fetch them.


The bowtie2 alignment is modified to output as well the unmapped reads using `--un-conc-gz` flag. The unmapped reads and total reads are then aligned to the reference, mitochondrial and chloroplast genome. The output of the mapping is then analysed with `samtools flagstat` and hopefully nicely ploted thanks to **MultiQc**.


# Content of the repository

- **Snakefile** containing the targeted output and the rules to generate them from the input files.

- **config/** , folder containing the configuration files making the Snakefile adaptable to any input files, genome and parameter for the rules. Adapt the config file and its reference in the Snakefile.

- **Fastq/**, folder containing subsetted paired-end fastq files used to test locally the pipeline. Generated using [Seqtk](https://github.com/lh3/seqtk): `seqtk sample -s100 read1.fq 5000 > sub1.fqseqtk sample -s100 read2.fq 5000 > sub2.fq`. RAW fastq or fastq.gz files should be placed here before running the pipeline.

- **envs/**, folder containing the environment needed for the Snakefile to run. To use Snakemake, it is required to create and activate an environment containing snakemake (here : envs/global_env.yaml )

- **units.tsv**, is a tab separated value files containing information about the experiment name, the condition of the experiment (control or treatment) and the path to the fastq files relative to the **Snakefile**. **Change this file according to your samples.**


# Usage

## Conda environment

First, you need to create an environment for the use of Snakemake with [Conda package manager](https://conda.io/docs/using/envs.html).
1. Create a virtual environment named "chipseq" from the `global_env.yaml` file with the following command: `conda env create --name chipseq --file ~/envs/global_env.yaml`
2. Then, activate this virtual environment with `source activate chipseq`

The Snakefile will then take care of installing and loading the packages and softwares required by each step of the pipeline.

## Configuration file
The `~/configs/config_tomato_sub.yaml` file specifies the sample list, the genomic reference fasta file to use, the directories to use, etc. This file is then used to build parameters in the main `Snakefile`.

## Snakemake execution
The Snakemake pipeline/workflow management system reads a master file (often called `Snakefile`) to list the steps to be executed and defining their order.
It has many rich features. Read more [here](https://snakemake.readthedocs.io/en/stable/).

## Samples
Samples are listed in the `units.tsv` file and will be used by the Snakefile automatically. Change the name, the conditions accordingly.

## Dry run
Use the command `snakemake -np` to perform a dry run that prints out the rules and commands.

## Real run
Simply type `Snakemake --use-conda` and provide the number of cores with `--cores 10` for ten cores for instance.
For cluster execution, please refer to the [Snakemake reference](https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution).
Please pay attention to `--use-conda`, it is required for the installation and loading of the dependencies used by the rules of the pipeline.
To run the pipeline, from the folder containing the Snakefile run the

# Main outputs
The main output are for now  **fastqc**, **bed** and **bigwig files**. Optionals outputs of the pipelines are **bamCompare**, **bedgraph** and **bed files for broad peaks calling**.
With the implementation of deeptools rules, the pipeline produces as well those outputs:

- **plotFingerprint** contains interesting figures that answer the question: **"Did my ChIP work???"** . Explanation of the plot and the options available can be found [here](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html)

- **PLOTCORRELATION** folder contain pdf files displaying the correlation between the samples tested in the ChIP experiment, many options in the plotcorrelation rules can be changed via the configuration file. More information about this plot can be found [here](https://deeptools.readthedocs.io/en/develop/content/tools/plotCorrelation.html)

- **HEATMAP** folder contain pdf files displaying the content of the matrix produced by the `computeMatrix` rule under the form of a heatmap. Many option for the `computeMatrix` and the `plotHeatmap` rules can be changed in the configuration file. More information about this figure can be found [here](https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html).

- **plotProfile** folder contain pdf files displaying profile plot for scores over sets of genomic region, again the genomic region are define in the matrix made previously. Again there are many options to change the plot and more information can be found [here](https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html)
