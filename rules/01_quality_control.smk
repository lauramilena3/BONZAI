rule countReads_gz:
    input:
        fastq=expand(
            dirs_dict["RAW_DATA_DIR"] + "/{sample}_L00{lane}_R{read}_001.fastq.gz",
            sample=SAMPLES, lane=["1", "2", "3", "4"], read=READ_TYPES
        ),
    output:
        counts=expand(
            dirs_dict["RAW_DATA_DIR"] + "/{sample}_L00{lane}_R{read}_read_count.txt",
            sample=SAMPLES, lane=["1", "2", "3", "4"], read=READ_TYPES
        ),
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
        raw_fastq=dirs_dict["RAW_DATA_DIR"] + "/{sample}_L00{lane}_R{read}_001.fastq.gz"
    output:
        html=temp(dirs_dict["RAW_DATA_DIR"] + "/{sample}_L00{lane}_R{read}_fastqc.html"),
        zipped=dirs_dict["RAW_DATA_DIR"] + "/{sample}_L00{lane}_R{read}_fastqc.zip"
    message:
        "Performing FastQC on raw data: {input.raw_fastq}"
    conda:
        dirs_dict["ENVS_DIR"] + "/QC.yaml"
    benchmark:
        dirs_dict["BENCHMARKS"] + "/01_QC/{sample}_L00{lane}_R{read}_pre_qc.tsv"
    shell:
        """
        fastqc {input.raw_fastq} --outdir {dirs_dict["RAW_DATA_DIR"]}
        """

rule fastQC_post:
    input:
        raw_fastq=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_L00{lane}_R{read}_001.fastq.gz"
    output:
        html=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_L00{lane}_R{read}_fastqc.html"),
        zipped=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_L00{lane}_R{read}_fastqc.zip"
    message:
        "Performing FastQC on cleaned data: {input.raw_fastq}"
    conda:
        dirs_dict["ENVS_DIR"] + "/QC.yaml"
    benchmark:
        dirs_dict["BENCHMARKS"] + "/01_QC/{sample}_L00{lane}_R{read}_post_qc.tsv"
    shell:
        """
        fastqc {input.raw_fastq} --outdir {dirs_dict["CLEAN_DATA_DIR"]}
        """

rule preMultiQC:
    input:
        zipped=expand(dirs_dict["RAW_DATA_DIR"] + "/{sample}_L00{lane}_R{read}_fastqc.zip",
                      sample=SAMPLES, lane=["1", "2", "3", "4"], read=READ_TYPES)
    output:
        multiqc=dirs_dict["QC_DIR"] + "/preQC_illumina_report.html",
        multiqc_txt=dirs_dict["QC_DIR"] + "/preQC_illumina_report_data/multiqc_fastqc.txt"
    params:
        fastqc_dir=dirs_dict["RAW_DATA_DIR"],
        html_name="preQC_illumina_report.html",
        multiqc_dir=dirs_dict["QC_DIR"]
    message:
        "Generating MultiQC report (Pre-QC)"
    conda:
        dirs_dict["ENVS_DIR"] + "/QC.yaml"
    benchmark:
        dirs_dict["BENCHMARKS"] + "/01_QC/multiqc_pre.tsv"
    shell:
        """
        multiqc -f {params.fastqc_dir} -o {params.multiqc_dir} -n {params.html_name}
        """

rule postMultiQC:
    input:
        zipped_forward=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_L00{lane}_R1_fastqc.zip",
                              sample=SAMPLES, lane=["1", "2", "3", "4"]),
        zipped_reverse=expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_L00{lane}_R2_fastqc.zip",
                              sample=SAMPLES, lane=["1", "2", "3", "4"])
    output:
        multiqc=dirs_dict["QC_DIR"] + "/postQC_illumina_report.html",
        multiqc_txt=dirs_dict["QC_DIR"] + "/postQC_illumina_report_data/multiqc_fastqc.txt"
    params:
        fastqc_dir=dirs_dict["CLEAN_DATA_DIR"],
        html_name="postQC_illumina_report.html",
        multiqc_dir=dirs_dict["QC_DIR"]
    message:
        "Generating MultiQC report (Post-QC)"
    conda:
        dirs_dict["ENVS_DIR"] + "/QC.yaml"
    benchmark:
        dirs_dict["BENCHMARKS"] + "/01_QC/multiqc_post.tsv"
    shell:
        """
        multiqc -f {params.fastqc_dir}/*zip -o {params.multiqc_dir} -n {params.html_name}
        """