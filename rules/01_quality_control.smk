rule countReads_gz:
    input:
        fastq=expand(dirs_dict["RAW_DATA_DIR"] + "/{sample}_L00{lane}_{read}_001.fastq.gz", sample=SAMPLES, lane=["1", "2", "3", "4"], read=READ_TYPES),
    output:
        counts=expand(dirs_dict["RAW_DATA_DIR"] + "/{sample}_L00{lane}_{read}_read_count.txt", sample=SAMPLES, lane=["1", "2", "3", "4"], read=READ_TYPES),
    message:
        "Counting reads in fastq.gz file",
    conda:
        dirs_dict["ENVS_DIR"] + "/QC.yaml",
    shell:
        """
        echo $(( $(zgrep -Ec "$" {input.fastq}) / 4 )) > {output.counts}
        """

rule fastQC_pre:
    input:
        raw_fastq=expand(dirs_dict["RAW_DATA_DIR"] + "/{sample}_L00{lane}_{read}_001.fastq.gz", lane=["1", "2", "3", "4"], read=READ_TYPES),
    output:
        html=temp(dirs_dict["RAW_DATA_DIR"] + "/{sample}_{read}_fastqc.html"),
        zipped=(dirs_dict["RAW_DATA_DIR"] + "/{sample}_{read}_fastqc.zip"),
    message:
        "Performing fastQC statistics on raw data",
    conda:
        dirs_dict["ENVS_DIR"] + "/QC.yaml",
    benchmark:
        dirs_dict["BENCHMARKS"] + "/01_QC/{sample}_pre_qc.tsv",
    shell:
        """
        fastqc {input.raw_fastq}
        """

rule fastQC_post:
    input:
        raw_fastq=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_L00{lane}_{read}_001.fastq.gz", lane=["1", "2", "3", "4"], read=READ_TYPES),
    output:
        html=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_{read}_fastqc.html"),
        zipped=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_{read}_fastqc.zip"),
    message:
        "Performing fastQC statistics on cleaned data",
    conda:
        dirs_dict["ENVS_DIR"] + "/QC.yaml",
    benchmark:
        dirs_dict["BENCHMARKS"] + "/01_QC/{sample}_post_qc.tsv",
    shell:
        """
        fastqc {input.raw_fastq}
        """

rule preMultiQC:
    input:
        zipped=expand(dirs_dict["RAW_DATA_DIR"] + "/{sample}_{read}_fastqc.zip", sample=SAMPLES, read=READ_TYPES),
    output:
        multiqc=dirs_dict["QC_DIR"] + "/preQC_illumina_report.html",
        multiqc_txt=dirs_dict["QC_DIR"] + "/preQC_illumina_report_data/multiqc_fastqc.txt",
    params:
        fastqc_dir=dirs_dict["RAW_DATA_DIR"],
        html_name="preQC_illumina_report.html",
        multiqc_dir=dirs_dict["QC_DIR"],
    message:
        "Generating MultiQC report (Pre-QC)",
    conda:
        dirs_dict["ENVS_DIR"] + "/QC.yaml",
    benchmark:
        dirs_dict["BENCHMARKS"] + "/01_QC/multiqc_pre.tsv",
    shell:
        """
        multiqc -f {params.fastqc_dir} -o {params.multiqc_dir} -n {params.html_name}
        """

rule postMultiQC:
    input:
        zipped_forward=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean_fastqc.zip", sample=SAMPLES),
        zipped_reverse=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean_fastqc.zip", sample=SAMPLES),
        zipped_unpaired=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean_fastqc.zip", sample=SAMPLES),
    output:
        multiqc=dirs_dict["QC_DIR"] + "/postQC_illumina_report.html",
        multiqc_txt=dirs_dict["QC_DIR"] + "/postQC_illumina_report_data/multiqc_fastqc.txt",
    params:
        fastqc_dir=dirs_dict["CLEAN_DATA_DIR"],
        html_name="postQC_illumina_report.html",
        multiqc_dir=dirs_dict["QC_DIR"],
    message:
        "Generating MultiQC report (Post-QC)",
    conda:
        dirs_dict["ENVS_DIR"] + "/QC.yaml",
    benchmark:
        dirs_dict["BENCHMARKS"] + "/01_QC/multiqc_post.tsv",
    shell:
        """
        multiqc -f {params.fastqc_dir}/*zip -o {params.multiqc_dir} -n {params.html_name}
        """