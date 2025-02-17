rule countReads_gz:
    input:
        fastq=config["dirs"]["RAW_DATA_DIR"] + "/{sample}_L00{lane}_{read}_001.fastq.gz"
    output:
        counts=config["dirs"]["RAW_DATA_DIR"] + "/{sample}_L00{lane}_{read}_read_count.txt"
    log:
        count_log=config["dirs"]["BENCHMARKS"] + "/01_QC/{sample}_L00{lane}_{read}_read_count.log"
    message:
        "Counting reads in fastq.gz file: {input.fastq}"
    conda:
        config["dirs"]["ENVS_DIR"] + "/QC.yaml"
    benchmark:
        config["dirs"]["BENCHMARKS"] + "/01_QC/{sample}_L00{lane}_{read}_read_count.tsv"
    shell:
        """
        echo "Checking file: {input.fastq}" > {log.count_log}
        ls -l {input.fastq} >> {log.count_log}
        zgrep -c "^" {input.fastq} | awk '{{print $1/4}}' > {output.counts} 2>> {log.count_log}
        """

rule fastQC_pre:
    input:
        raw_fastq=config["dirs"]["RAW_DATA_DIR"] + "/{sample}_L00{lane}_{read}_001.fastq.gz"
    output:
        html=temp(config["dirs"]["RAW_DATA_DIR"] + "/{sample}_L00{lane}_{read}_fastqc.html"),
        zipped=config["dirs"]["RAW_DATA_DIR"] + "/{sample}_L00{lane}_{read}_fastqc.zip"
    message:
        "Performing FastQC on raw data: {input.raw_fastq}"
    conda:
        config["dirs"]["ENVS_DIR"] + "/QC.yaml"
    shell:
        """
        fastqc {input.raw_fastq} --outdir {config["dirs"]["RAW_DATA_DIR"]}
        """

rule fastQC_post:
    input:
        raw_fastq=config["dirs"]["QC_DIR"] + "/{sample}_L00{lane}_{read}_001.fastq.gz"
    output:
        html=temp(config["dirs"]["QC_DIR"] + "/{sample}_L00{lane}_{read}_fastqc.html"),
        zipped=config["dirs"]["QC_DIR"] + "/{sample}_L00{lane}_{read}_fastqc.zip"
    message:
        "Performing FastQC on cleaned data: {input.raw_fastq}"
    conda:
        config["dirs"]["ENVS_DIR"] + "/QC.yaml"
    shell:
        """
        fastqc {input.raw_fastq} --outdir {config["dirs"]["QC_DIR"]}
        """