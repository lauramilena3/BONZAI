rule countReads_gz:
    input:
        fastq=dirs_dict["RAW_DATA_DIR"] + "/{sample}_L00{lane}_{read}_001.fastq.gz"
    output:
        counts=dirs_dict["RAW_DATA_DIR"] + "/{sample}_L00{lane}_{read}_read_count.txt"
    log:
        count_log=dirs_dict["BENCHMARKS"] + "/01_QC/{sample}_L00{lane}_{read}_read_count.log"
    message:
        "Counting reads in fastq.gz file: {input.fastq}"
    conda:
        dirs_dict["ENVS_DIR"] + "/QC.yaml"
    benchmark:
        dirs_dict["BENCHMARKS"] + "/01_QC/{sample}_L00{lane}_{read}_read_count.tsv"
    shell:
        """
        echo "Checking file: {input.fastq}" > {log.count_log}
        ls -l {input.fastq} >> {log.count_log}
        zgrep -c "^" {input.fastq} | awk '{{print $1/4}}' > {output.counts} 2>> {log.count_log}
        """

rule fastQC_pre:
    input:
        raw_fastq=dirs_dict["RAW_DATA_DIR"] + "/{sample}_L00{lane}_{read}_001.fastq.gz"
    output:
        html=temp(dirs_dict["RAW_DATA_DIR"] + "/{sample}_L00{lane}_{read}_fastqc.html"),
        zipped=dirs_dict["RAW_DATA_DIR"] + "/{sample}_L00{lane}_{read}_fastqc.zip"
    log:
        fastqc_log=dirs_dict["BENCHMARKS"] + "/01_QC/{sample}_L00{lane}_{read}_fastqc.log"
    benchmark:
        dirs_dict["BENCHMARKS"] + "/01_QC/{sample}_L00{lane}_{read}_pre_qc.tsv"
    message:
        "Performing FastQC on raw data: {input.raw_fastq}"
    conda:
        dirs_dict["ENVS_DIR"] + "/QC.yaml"
    shell:
        """
        fastqc {input.raw_fastq} --outdir {dirs_dict["RAW_DATA_DIR"]} &> {log.fastqc_log}
        """

rule fastQC_post:
    input:
        raw_fastq=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_L00{lane}_{read}_001.fastq.gz"
    output:
        html=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_L00{lane}_{read}_fastqc.html"),
        zipped=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_L00{lane}_{read}_fastqc.zip"
    log:
        fastqc_log=dirs_dict["BENCHMARKS"] + "/01_QC/{sample}_L00{lane}_{read}_fastqc.log"
    benchmark:
        dirs_dict["BENCHMARKS"] + "/01_QC/{sample}_L00{lane}_{read}_post_qc.tsv"
    message:
        "Performing FastQC on cleaned data: {input.raw_fastq}"
    conda:
        dirs_dict["ENVS_DIR"] + "/QC.yaml"
    shell:
        """
        fastqc {input.raw_fastq} --outdir {dirs_dict["CLEAN_DATA_DIR"]} &> {log.fastqc_log}
        """

rule preMultiQC:
    input:
        zipped=expand(dirs_dict["RAW_DATA_DIR"] + "/{sample}_L00{lane}_{read}_fastqc.zip",
                      sample=SAMPLES, lane=["1", "2", "3", "4"], read=["R1", "R2"])
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
