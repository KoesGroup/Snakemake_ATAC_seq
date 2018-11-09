rule trimmomatic:
    input:
        reads = get_fastq,
        adapters = config["adapters"]
    output:
        forward_reads   = WORKING_DIR + "trimmed/{sample}_forward.fastq.gz",
        reverse_reads   = WORKING_DIR + "trimmed/{sample}_reverse.fastq.gz",
        forwardUnpaired = temp(WORKING_DIR + "trimmed/{sample}_forward_unpaired.fastq.gz"),
        reverseUnpaired = temp(WORKING_DIR + "trimmed/{sample}_reverse_unpaired.fastq.gz")
    message: "trimming {wildcards.sample} reads"
    log:
        RESULT_DIR + "logs/trimmomatic/{sample}.log"
    params :
        seedMisMatches =            str(config['trimmomatic']['seedMisMatches']),
        palindromeClipTreshold =    str(config['trimmomatic']['palindromeClipTreshold']),
        simpleClipThreshhold =      str(config['trimmomatic']['simpleClipThreshold']),
        LeadMinTrimQual =           str(config['trimmomatic']['LeadMinTrimQual']),
        TrailMinTrimQual =          str(config['trimmomatic']['TrailMinTrimQual']),
        windowSize =                str(config['trimmomatic']['windowSize']),
        avgMinQual =                str(config['trimmomatic']['avgMinQual']),
        minReadLen =                str(config['trimmomatic']['minReadLength']),
        phred = 		            str(config["trimmomatic"]["phred"])
    threads: 10
    conda:
        "../envs/trimmomatic_env.yaml"
    shell:
        "trimmomatic PE {params.phred} -threads {threads} "
        "{input.reads} "
        "{output.forward_reads} "
        "{output.forwardUnpaired} "
        "{output.reverse_reads} "
        "{output.reverseUnpaired} "
        "ILLUMINACLIP:{input.adapters}:{params.seedMisMatches}:{params.palindromeClipTreshold}:{params.simpleClipThreshhold} "
        "LEADING:{params.LeadMinTrimQual} "
        "TRAILING:{params.TrailMinTrimQual} "
        "SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} "
        "MINLEN:{params.minReadLen} &>{log}"

rule fastqc:
    input:
        fwd = WORKING_DIR + "trimmed/{sample}_forward.fastq.gz",
        rev = WORKING_DIR + "trimmed/{sample}_reverse.fastq.gz"
    output:
        fwd = RESULT_DIR + "fastqc/{sample}_forward_fastqc.zip",
        rev = RESULT_DIR + "fastqc/{sample}_reverse_fastqc.zip"
    log:
        RESULT_DIR + "logs/fastqc/{sample}.fastqc.log"
    params:
        RESULT_DIR + "fastqc/"
    message:
        "---Quality check of trimmed {wildcards.sample} sample with FASTQC"
    conda:
        "../envs/fastqc_env.yaml"
    shell:
        "fastqc --outdir={params} {input.fwd} {input.rev} &>{log}"

rule index:
    input:
        WORKING_DIR + "genome.fasta"
    output:
        [WORKING_DIR + "genome." + str(i) + ".bt2" for i in range(1,5)],
        WORKING_DIR + "genome.rev.1.bt2",
        WORKING_DIR + "genome.rev.2.bt2"
    message:"indexing genome"
    params:
        WORKING_DIR + "genome"
    threads: 10
    conda:
        "../envs/samtools_bowtie_env.yaml"
    shell:"bowtie2-build --threads {threads} {input} {params}"

rule index_chloro:
    input:
        "SL_chloroplaste_sequence.fasta"
    output:
        [WORKING_DIR + "chloroplast." + str(i) + ".bt2" for i in range(1,5)],
        WORKING_DIR + "chloroplast.rev.1.bt2",
        WORKING_DIR + "chloroplast.rev.2.bt2"
    message:"indexing chloroplast genome"
    params:
        WORKING_DIR + "chloroplast"
    threads: 10
    conda:
        "../envs/samtools_bowtie_env.yaml"
    shell:"bowtie2-build --threads {threads} {input} {params}"

rule index_mito:
    input:
        "S_lycopersicum_mitochondrion_v1.5_genomic.fna"
    output:
        [WORKING_DIR + "mitochondrialgenome." + str(i) + ".bt2" for i in range(1,5)],
        WORKING_DIR + "mitochondrialgenome.rev.1.bt2",
        WORKING_DIR + "mitochondrialgenome.rev.2.bt2"
    message:"indexing mitochondrial genome"
    params:
        WORKING_DIR + "mitochondrialgenome"
    threads: 10
    conda:
        "../envs/samtools_bowtie_env.yaml"
    shell:"bowtie2-build --threads {threads} {input} {params}"

rule align:
    input:
        forward         = WORKING_DIR + "trimmed/{sample}_forward.fastq.gz",
        reverse         = WORKING_DIR + "trimmed/{sample}_reverse.fastq.gz",
        forwardUnpaired = WORKING_DIR + "trimmed/{sample}_forward_unpaired.fastq.gz",
        reverseUnpaired = WORKING_DIR + "trimmed/{sample}_reverse_unpaired.fastq.gz",
        index           = [WORKING_DIR + "genome." + str(i) + ".bt2" for i in range(1,5)]
    output:
        mapped          = WORKING_DIR + "mapped/{sample}.bam",
        unmapped        = [WORKING_DIR + "unmapped/{sample}.fq." + str(i) +".gz" for i in range(1,2)],
        bai             = WORKING_DIR + "mapped/{sample}.sorted.bam.bai",
        sorted          = WORKING_DIR + "mapped/{sample}.sorted.bam"
    message: "Mapping files {wildcards.sample}"
    params:
        bowtie          = " ".join(config["bowtie2"]["params"].values()), #take argument separated as a list separated with a space
        index           = WORKING_DIR + "genome",
        unmapped        = WORKING_DIR + "unmapped/{sample}.fq.gz"
    threads: 10
    conda:
        "../envs/samtools_bowtie_env.yaml"
    log:
        RESULT_DIR + "logs/bowtie/{sample}.log"
    shell:
        """
        bowtie2 {params.bowtie} --threads {threads} -x {params.index} -1 {input.forward} -2 {input.reverse} -U {input.forwardUnpaired},{input.reverseUnpaired} --un-conc-gz {params.unmapped} | samtools view -Sb - > {output.mapped} 2>{log}
        samtools sort -o {output.sorted} {output.mapped} &>{log}
        samtools index {output.sorted}
        """    
rule align_chloro:
    input:
        forward         = WORKING_DIR + "trimmed/{sample}_forward.fastq.gz",
        reverse         = WORKING_DIR + "trimmed/{sample}_reverse.fastq.gz",
        forwardUnpaired = WORKING_DIR + "trimmed/{sample}_forward_unpaired.fastq.gz",
        reverseUnpaired = WORKING_DIR + "trimmed/{sample}_reverse_unpaired.fastq.gz",
        index           = [WORKING_DIR + "chloroplast." + str(i) + ".bt2" for i in range(1,5)],
        nonmapped_for   = WORKING_DIR + "unmapped/{sample}.fq.1.gz",
        nonmapped_rev   = WORKING_DIR + "unmapped/{sample}.fq.2.gz"
    output:
        mapped          = WORKING_DIR + "mapped/chloroplast_mapped/{sample}.bam",
        unmapped        = WORKING_DIR + "mapped/chloroplast_unmapped/{sample}.bam"
    conda:
        "../envs/samtools_bowtie_env.yaml"
    message:
        "Mapping files {wildcards.sample} to chloroplast DNA"
    params:
        bowtie          = " ".join(config["bowtie2"]["params"].values()), #take argument separated as a list separated with a space
        index           = WORKING_DIR + "chloroplast"
    threads: 10
    shell:
        """
        bowtie2 {params.bowtie} --threads {threads} -x {params.index} -1 {input.forward} -2 {input.reverse} -U {input.forwardUnpaired},{input.reverseUnpaired} | samtools view -Sb - > {output.mapped}
        bowtie2 {params.bowtie} --threads {threads} -x {params.index} -1 {input.nonmapped_for} -2 {input.nonmapped_rev}| samtools view -Sb - > {output.unmapped}
        """

rule align_mito:
    input:
        forward         = WORKING_DIR + "trimmed/{sample}_forward.fastq.gz",
        reverse         = WORKING_DIR + "trimmed/{sample}_reverse.fastq.gz",
        forwardUnpaired = WORKING_DIR + "trimmed/{sample}_forward_unpaired.fastq.gz",
        reverseUnpaired = WORKING_DIR + "trimmed/{sample}_reverse_unpaired.fastq.gz",
        index           = [WORKING_DIR + "mitochondrialgenome." + str(i) + ".bt2" for i in range(1,5)],
        nonmapped_for   = WORKING_DIR + "unmapped/{sample}.fq.1.gz",
        nonmapped_rev   = WORKING_DIR + "unmapped/{sample}.fq.2.gz"
    output:
        mapped          = WORKING_DIR + "mapped/mitochondrial_mapped/{sample}.bam",
        unmapped        = WORKING_DIR + "mapped/mitochondrial_unmapped/{sample}.bam"
    conda:
        "../envs/samtools_bowtie_env.yaml"
    message:
        "Mapping files {wildcards.sample} to mitochondrial genome"
    params:
        bowtie          = " ".join(config["bowtie2"]["params"].values()), #take argument separated as a list separated with a space
        index           = WORKING_DIR + "mitochondrialgenome"
    threads: 10
    shell:
        """
        bowtie2 {params.bowtie} --threads {threads} -x {params.index} -1 {input.forward} -2 {input.reverse} -U {input.forwardUnpaired},{input.reverseUnpaired} | samtools view -Sb - > {output.mapped}
        bowtie2 {params.bowtie} --threads {threads} -x {params.index} -1 {input.nonmapped_for} -2 {input.nonmapped_rev}| samtools view -Sb - > {output.unmapped}
        """

rule flagstat_genome:
    input:
        WORKING_DIR + "mapped/{sample}.bam"
    output:
        RESULT_DIR + "logs/flagstat/genome/{sample}.bam.flagstat"
    message:
        "Analyzing mapping information of {wildcards.sample}"
    params:
        jobname = "{sample}"
    conda:
        "../envs/samtools_bowtie_env.yaml"
    shell:
        "samtools flagstat {input} > {output}"

rule flagstat_mitochondrial:
    input:
        WORKING_DIR + "mapped/mitochondrial_mapped/{sample}.bam"
    output:
        RESULT_DIR + "logs/flagstat/mitochondrial/{sample}.bam.flagstat"
    message:
        "Analyzing mapping information of {wildcards.sample}"
    params:
        jobname = "{sample}"
    conda:
        "../envs/samtools_bowtie_env.yaml"
    shell:
        "samtools flagstat {input} > {output}"

rule flagstat_chloroplast:
    input:
        WORKING_DIR + "mapped/chloroplast_mapped/{sample}.bam"
    output:
        RESULT_DIR + "logs/flagstat/chloroplast/{sample}.bam.flagstat"
    message:
        "Analyzing mapping information of {wildcards.sample}"
    params:
        jobname = "{sample}"
    conda:
        "../envs/samtools_bowtie_env.yaml"
    shell:
        "samtools flagstat {input} > {output}"

rule alignmentsieve:
    input:
        bam = WORKING_DIR + "mapped/{sample}.sorted.bam",
        bai = WORKING_DIR + "mapped/{sample}.sorted.bam.bai"
    output:
        bam              = RESULT_DIR + "mapped/{sample}.shifted.rmdup.bam",
        filteredOutReads = RESULT_DIR + "mapped/filteredOutReads/{sample}.bam",
        filterMetrics    = RESULT_DIR + "logs/alignmentsieve/{sample}.txt" 
    log:
        RESULT_DIR + "logs/alignmentsieve/{sample}.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "alignmentSieve --bam {input.bam} --outFile {output.bam} --ATACshift --filteredOutReads {output.filteredOutReads} --filterMetrics {output.filterMetrics} --ignoreDuplicates 2>{log}"

rule sort:
    input:
        RESULT_DIR + "mapped/{sample}.shifted.rmdup.bam"
    output:
        bam = RESULT_DIR + "mapped/{sample}.shifted.rmdup.sorted.bam",
        bai = RESULT_DIR + "mapped/{sample}.shifted.rmdup.sorted.bam.bai"
    message:"sorting {wildcards.sample} bam file"
    threads: 10
    log:
        RESULT_DIR + "logs/samtools/{sample}.sort.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input} &>{log}
        samtools index {output.bam}
        """