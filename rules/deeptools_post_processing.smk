rule bamCoverage:
    input:
        RESULT_DIR + "mapped/{sample}.shifted.rmdup.sorted.bam"
    output:
        RESULT_DIR + "bigwig/{sample}.bw"
    message:
        "Converting {wildcards.sample} bam into bigwig file"
    log:
        RESULT_DIR + "logs/deeptools/{sample}_bamtobigwig.log"
    params:
        EFFECTIVEGENOMESIZE     = str(config["bamCoverage"]["params"]["EFFECTIVEGENOMESIZE"]), #take argument separated as a list separated with a space
        EXTENDREADS             = str(config["bamCoverage"]["params"]["EXTENDREADS"]),
        binSize                 = str(config['bamCoverage']["params"]['binSize']),
        normalizeUsing          = str(config['bamCoverage']["params"]['normalizeUsing']),
        ignoreForNormalization  = str(config['bamCoverage']["params"]['ignoreForNormalization']),
        smoothLength            = str(config['bamCoverage']["params"]['smoothLength'])
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bamCoverage --bam {input} \
        -o {output} \
        --effectiveGenomeSize {params.EFFECTIVEGENOMESIZE} \
        --extendReads {params.EXTENDREADS} \
        --binSize {params.binSize} \
        --smoothLength {params.smoothLength} \
        --ignoreForNormalization {params.ignoreForNormalization} \
        &>{log}"

rule bamcompare:
    input:
        treatment   = RESULT_DIR + "mapped/{treatment}.shifted.rmdup.sorted.bam",              #input requires an indexed bam file
        control     = RESULT_DIR + "mapped/{control}.shifted.rmdup.sorted.bam"                   #input requires an indexed bam file
    output:
        bigwig = RESULT_DIR + "bamcompare/log2_{treatment}_{control}.bamcompare.bw"
    message:
        "Running bamCompare for {wildcards.treatment} and {wildcards.control}"
    log:
        RESULT_DIR + "logs/deeptools/log2_{treatment}_{control}.bamcompare.bw.log"
    conda:
        "../envs/deeptools.yaml"
    params:
        binSize             = str(config['bamcompare']['binSize']),
        normalizeUsing      = str(config['bamcompare']['normalizeUsing']),
        EFFECTIVEGENOMESIZE = str(config["bamcompare"]["EFFECTIVEGENOMESIZE"]),
        operation           = str(config['bamcompare']['operation']),
        smoothLength        = str(config['bamcompare']['smoothLength']),
        ignoreForNormalization = str(config['bamcompare']['ignoreForNormalization']),
        scaleFactorsMethod  = str(config['bamcompare']['scaleFactorsMethod'])
    shell:
        "bamCompare -b1 {input.treatment} \
        -b2 {input.control}  \
        --binSize {params.binSize} \
        -o {output.bigwig} \
        --normalizeUsing {params.normalizeUsing} \
        --operation {params.operation} \
        --smoothLength {params.smoothLength} \
        --ignoreForNormalization {params.ignoreForNormalization} \
        --scaleFactorsMethod {params.scaleFactorsMethod} \
        &>{log}"

rule multiBamSummary:
    input:
        lambda wildcards: expand(RESULT_DIR + "mapped/{sample}.shifted.rmdup.sorted.bam", sample = SAMPLES)
    output:
        RESULT_DIR + "multiBamSummary/MATRIX.npz"
    message:
        "Computing the read coverage into a numpy array "
    threads: 10
    params:
        binSize     = str(config['multiBamSummary']['binSize'])
    log:
        RESULT_DIR + "logs/deeptools/multibamsummary/MATRIX.log"

    shell:
        "multiBamSummary bins \
        --bamfiles {input} \
        --numberOfProcessors {threads}\
        --binSize {params.binSize} \
        --centerReads \
        --extendReads \
        -o {output} \
        2> {log}"


rule plotCorrelation:
    input:
        RESULT_DIR + "multiBamSummary/MATRIX.npz"
    output:
        RESULT_DIR + "plotCorrelation/MATRIX.png"
    log:
        RESULT_DIR + "logs/deeptools/plotcorrelation/MATRIX.log"
    params:
        corMethod  = str(config['plotCorrelation']['corMethod']),
        whatToPlot = str(config['plotCorrelation']['whatToPlot']),
        color      = str(config['plotCorrelation']['color'])
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotCorrelation \
                    --corData {input} \
                    --corMethod {params.corMethod} \
                    --whatToPlot {params.whatToPlot} \
                    --skipZeros \
                    --colorMap {params.color} \
                    --plotFile {output} \
                    --plotNumbers \
                    2> {log}"

rule computeMatrix:
    input:
        bigwig = RESULT_DIR + "bigwig/{sample}.bw",
        bed    = WORKING_DIR + "gene_model.gtf"
    output:
        RESULT_DIR + "computematrix/{sample}.TSS.gz"
    threads: 10
    params:
        binSize = str(config['computeMatrix']['binSize'])
    conda:
        "../envs/deeptools.yaml"
    log:
        RESULT_DIR + "logs/deeptools/computematrix/{sample}.log"
    shell:
        "computeMatrix \
        reference-point \
        --referencePoint TSS \
        -S {input.bigwig} \
        -R {input.bed} \
        --afterRegionStartLength 3000 \
        --beforeRegionStartLength 3000 \
        --numberOfProcessors {threads} \
        --binSize {params.binSize} \
        -o {output} \
        2> {log}"

rule plotHeatmap:
    input:
        RESULT_DIR + "computematrix/{sample}.TSS.gz"
    output:
        RESULT_DIR + "heatmap/{sample}.pdf"
    params:
        kmeans = str(config['plotHeatmap']['kmeans']),
        color  = str(config['plotHeatmap']['color']),
        plot   = str(config['plotHeatmap']['plot']),
        cluster = RESULT_DIR + "heatmap/{sample}.bed"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotHeatmap \
        --matrixFile {input} \
        --outFileName {output} \
        --kmeans {params.kmeans} \
        --colorMap {params.color} \
        --legendLocation best \
        --outFileSortedRegions {params.cluster}"

rule plotFingerprint:
    input:
        expand(RESULT_DIR + "mapped/{sample}.shifted.rmdup.sorted.bam", sample = SAMPLES)
    output:
        pdf = RESULT_DIR + "plotFingerprint/Fingerplot.pdf"
    params:
        EXTENDREADS  = str(config["bamCoverage"]["params"]["EXTENDREADS"]),
        binSize      = str(config['bamCoverage']["params"]['binSize'])
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotFingerprint \
        -b {input}\
        --extendReads {params.EXTENDREADS} \
        --binSize {params.binSize} \
        --plotFile {output}"

rule plotProfile:
    input:
        RESULT_DIR + "computematrix/{sample}.TSS.gz"
    output:
        pdf = RESULT_DIR + "plotProfile/{sample}.pdf",
        bed = RESULT_DIR + "plotProfile/{sample}.bed"
    params:
        kmeans      = str(config['plotProfile']['kmeans']),
        startLabel  = str(config['plotProfile']['startLabel']),
        endLabel    = str(config['plotProfile']['endLabel'])
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotProfile \
        --matrixFile {input} \
        --outFileName {output.pdf} \
        --outFileSortedRegions {output.bed} \
        --kmeans {params.kmeans} \
        --startLabel {params.startLabel} \
        --endLabel {params.endLabel}"

rule bamPEFragmentSize:
    input: expand(RESULT_DIR + "mapped/{sample}.shifted.rmdup.sorted.bam", sample = SAMPLES)
    output:
        png = RESULT_DIR + "bamPEFragmentSize/fragmentSize.png",
        RawFragmentLengths = RESULT_DIR + "bamPEFragmentSize/raw"
    conda:
        "../envs/deeptools.yaml"
    params:
        binSize = str(config['bamCoverage']["params"]['binSize']),
        #title   = "Fragment size of PE ATAC-seq data",
        maxFragmentLength = str(config['bamPEFragmentSize']['binSize'])   
    shell:
        "bamPEFragmentSize\
        --histogram {output.png} \
        --maxFragmentLength {params.maxFragmentLength} \
        --bamfiles {input} \
        --samplesLabel ATAC_1 ATAC_2 ATAC_3 ATAC_4 ATAC_5 ATAC_6 \
        --binSize {params.binSize} \
        --outRawFragmentLengths {output.RawFragmentLengths}"