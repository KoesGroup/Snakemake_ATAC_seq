# Snakemake file for ATAC-Seq PE analysis

###############
# Libraries
###############

import os
import pandas as pd
from snakemake.utils import validate, min_version
#############################################
# Configuration and sample sheets
#############################################

configfile: "config.yaml"

WORKING_DIR         = config["working_dir"]    # where you want to store your intermediate files (this directory will be cleaned up at the end)
RESULT_DIR          = config["result_dir"]      # what you want to keep

GENOME_FASTA_URL    = config["refs"]["genome_url"]
GENOME_FASTA_FILE   = os.path.basename(config["refs"]["genome_url"])
GFF_URL             = config["refs"]["gff_url"]
GFF_FILE            = os.path.basename(config["refs"]["gff_url"])
MT_GENOME_URL       = config["refs"]["mt_genome_url"]
MT_GENOME_FILE      = os.path.basename(config["refs"]["mt_genome_url"])

TOTALCORES          = 16                             #check this via 'grep -c processor /proc/cpuinfo'

###############
# Helper Functions
###############
def get_fastq(wildcards):
    return units.loc[(wildcards.sample), ["fq1", "fq2"]].dropna()

def get_samples_per_treatment(input_df="units.tsv",colsamples="sample",coltreatment="condition",treatment="control"):
    """This function returns a list of samples that correspond to the same experimental condition"""
    df = pd.read_table(input_df)
    df = df.loc[df[coltreatment] == treatment]
    filtered_samples = df[colsamples].tolist()
    return filtered_samples

def is_single_end(sample):
    """This function detect missing value in the column 2 of the units.tsv"""
    return pd.isnull(units.loc[(sample), "fq2"])

##############
# Samples and conditions
##############

units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)

SAMPLES = units.index.get_level_values('sample').unique().tolist()

CASES = get_samples_per_treatment(treatment="treatment")
CONTROLS = get_samples_per_treatment(treatment="control")

GROUPS = {
    "group1" : ["ATAC1", "ATAC2", "ATAC3"],
    "group2" : ["ATAC4", "ATAC5", "ATAC6"]
}                                           #I used this dictionnary to define the group of sample used in the multiBamSummary, might be improved a lot

##############
# Wildcards
##############
wildcard_constraints:
    sample = "[A-Za-z0-9]+"

wildcard_constraints:
    unit = "L[0-9]+"

##############
# Desired output
##############

FASTQC_REPORTS  =     expand(RESULT_DIR + "fastqc/{sample}_{pair}_fastqc.zip", sample=SAMPLES, pair={"forward", "reverse"})
BAM_INDEX       =     expand(RESULT_DIR + "mapped/{sample}.shifted.rmdup.sorted.bam.bai", sample=SAMPLES)
FLAGSTAT_GEN    =     expand(RESULT_DIR + "logs/flagstat/genome/{sample}.bam.flagstat", sample=SAMPLES)
FLAGSTAT_MITO   =     expand(RESULT_DIR + "logs/flagstat/mitochondrial/{sample}.bam.flagstat", sample=SAMPLES)
FLAGSTAT_CHLORO =     expand(RESULT_DIR + "logs/flagstat/chloroplast/{sample}.bam.flagstat", sample=SAMPLES)
BIGWIG          =     expand(RESULT_DIR + "bigwig/{sample}.bw", sample=SAMPLES)
BED_NARROW      =     expand(RESULT_DIR + "bed/{sample}_peaks.narrowPeak", sample=SAMPLES)
MULTIBAMSUMMARY =     RESULT_DIR + "multiBamSummary/MATRIX.npz"
PLOTCORRELATION =     RESULT_DIR + "plotCorrelation/MATRIX.png"
COMPUTEMATRIX   =     expand(RESULT_DIR + "computematrix/{sample}.{type}.gz", sample=SAMPLES, type={"TSS", "scale-regions"})
HEATMAP         =     expand(RESULT_DIR + "heatmap/{sample}.{type}.pdf", sample=SAMPLES, type={"TSS", "scale-regions"}) 
PLOTFINGERPRINT =     RESULT_DIR + "plotFingerprint/Fingerplot.pdf"
PLOTPROFILE_PDF =     expand(RESULT_DIR + "plotProfile/{sample}.{type}.pdf", sample=SAMPLES, type={"TSS", "scale-regions"})
PLOTPROFILE_BED =     expand(RESULT_DIR + "plotProfile/{sample}.{type}.bed", sample=SAMPLES, type={"TSS", "scale-regions"})
MULTIQC         =     "qc/multiqc.html"
FRAGMENTSIZE    =     RESULT_DIR + "bamPEFragmentSize/fragmentSize.png"
PLOTCOVERAGE    =     RESULT_DIR + "plotCoverage/Coverage.png"    

###############
# Final output
################
rule all:
    input:
        BAM_INDEX,
        FASTQC_REPORTS,
        BIGWIG,
        BED_NARROW,
        MULTIBAMSUMMARY,
        PLOTCORRELATION,
        COMPUTEMATRIX,
        HEATMAP,
        PLOTFINGERPRINT,
        PLOTPROFILE_PDF,
        PLOTPROFILE_BED,
        MULTIQC,
        FLAGSTAT_GEN,
        FLAGSTAT_MITO,
        FLAGSTAT_CHLORO,
        #FRAGMENTSIZE,
        PLOTCOVERAGE 
    message: "ATAC-seq pipeline succesfully run."		#finger crossed to see this message!

    shell:"#rm -rf {WORKING_DIR}"

###############
# Rules
###############

include : "rules/external_data.smk"
include : 'rules/pre_processing.smk'
include : "rules/macs2_peak_calling.smk"
include : "rules/deeptools_post_processing.smk"

# rule multiqc:
#     input:
#         RESULT_DIR + "logs/",
#         expand(RESULT_DIR + "fastqc/{sample}_{pair}_fastqc.zip", sample=SAMPLES, pair={"forward", "reverse"}),
#         RESULT_DIR + "plotCoverage/coverage.tab",
#         expand(RESULT_DIR + "bed/{sample}_peaks.xls", sample= SAMPLES)
#     output:
#         "multiqc_report.html/multiqc_report.html"
#     conda:
#         "envs/multiqc_env.yaml"
#     shell:
#         "multiqc {input} \
#         -o {output} \
#         -v \
#         -d \
#         -f "

         #RESULT_DIR + "bamPEFragmentSize/fragmentSize.tab",


rule multiqc:
    input:
        expand(RESULT_DIR + "fastqc/{sample}_{pair}_fastqc.zip", sample=SAMPLES, pair={"forward", "reverse"}),
        expand(RESULT_DIR + "bed/{sample}_peaks.xls", sample= SAMPLES)
    output:
        "qc/multiqc.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log"
    wrapper:
        "0.27.1/bio/multiqc"