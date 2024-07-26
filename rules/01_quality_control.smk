rule countReads_gz:
	input:
		fastq="{fastq_name}.fastq.gz",
	output:
		counts="{fastq_name}_read_count.txt",
	message:
		"Counting reads on fastq.gz file"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	shell:
		"""
		echo $(( $(zgrep -Ec "$" {input.fastq}) / 4 )) > {output.counts}
		"""

rule countReads:
	input:
		fastq="{fastq_name}.fastq",
	output:
		counts="{fastq_name}_read_count.txt",
	message:
		"Counting reads on fastq file"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	shell:
		"""
		echo $(( $(grep -Ec "$" {input.fastq}) / 4 )) > {output.counts}
		"""

rule fastQC_pre:
	input:
		raw_fastq=dirs_dict["RAW_DATA_DIR"] + "/{fastq_name}.fastq.gz"
	output:
		html=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{fastq_name}_fastqc.html"),
		zipped=(dirs_dict["CLEAN_DATA_DIR"] + "/{fastq_name}_fastqc.zip")
	message:
		"Performing fastqQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/qualityCheckIllumina/{fastq_name}_pre_qc.tsv"
	shell:
		"""
		fastqc {input}
		"""

rule fastQC_post:
	input:
		raw_fastq=dirs_dict["CLEAN_DATA_DIR"] + "/{fastq_name}.fastq.gz"
	output:
		html=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{fastq_name}_fastqc.html"),
		zipped=(dirs_dict["CLEAN_DATA_DIR"] + "/{fastq_name}_fastqc.zip")
	message:
		"Performing fastqQC statistics"
	conda:
		dirs_dict["ENVS_DIR"] + "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/qualityCheckIllumina/{fastq_name}_post_qc.tsv"
	shell:
		"""
		fastqc {input}
		"""

rule trim_adapters_quality_illumina_PE:
	input:
		forward_file=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['forward_tag']) + ".fastq.gz",
		reverse_file=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['reverse_tag']) + ".fastq.gz",
		adapters=dirs_dict["ADAPTERS_DIR"] + "/" + config['adapters_file']
	output:
		forward_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq.gz"),
		reverse_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq.gz"),
		forward_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_unpaired.fastq.gz"),
		reverse_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_unpaired.fastq.gz"),
		merged_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_merged_unpaired.fastq.gz"),
	message:
		"Trimming Illumina Adapters with Trimmomatic"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/trim_adapters_quality_illumina_PE/{sample}.tsv"
	threads: 8
	shell:
		"""
		trimmomatic PE -threads {threads} -phred33 {input.forward_file} {input.reverse_file} \
			{output.forward_paired} {output.forward_unpaired} {output.reverse_paired} {output.reverse_unpaired} \
			ILLUMINACLIP:{input.adapters}:2:30:10:1:true LEADING:{config[trimmomatic_leading]} TRAILING:{config[trimmomatic_trailing]} \
			SLIDINGWINDOW:{config[trimmomatic_window_size]}:{config[trimmomatic_window_quality]} MINLEN:{config[trimmomatic_minlen]}
		cat {output.forward_unpaired} {output.reverse_unpaired} > {output.merged_unpaired}
		"""

rule contaminants_KRAKEN:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq.gz"),
		merged_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_merged_unpaired.fastq.gz"),
		kraken_db=(config['kraken_db']),
		kraken_tools=(config['kraken_tools']),
	output:
		kraken_output_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output_paired.csv"),
		kraken_report_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired.csv"),
		kraken_output_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output_unpaired.csv"),
		kraken_report_unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_unpaired.csv"),
	params:
		kraken_db=config['kraken_db'],
	message:
		"Assesing contamination with kraken2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/kraken/{sample}_preliminary.tsv"
	threads: 8
	shell:
		"""
		kraken2 --db {params.kraken_db} --threads {threads} \
			--paired {input.forward_paired} {input.reverse_paired} \
			--output {output.kraken_output_paired} --report {output.kraken_report_paired}
		#UNPAIRED
		kraken2 --db {params.kraken_db} --threads {threads} {input.merged_unpaired}  \
			--output {output.kraken_output_unpaired} --report {output.kraken_report_unpaired}
		"""

rule remove_contaminants:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired.fastq.gz"),
		merged_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_merged_unpaired.fastq.gz"),
		kraken_output_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output_paired.csv"),
		kraken_output_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output_unpaired.csv"),
		kraken_report_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired.csv"),
		kraken_report_unpaired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_unpaired.csv"),
		kraken_tools=(config['kraken_tools']),
	output:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.fastq.gz"),
		unpaired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.fastq"),
	message:
		"Removing contaminants with Kraken"
	params:
		# unclassified_name_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken_paired_R#.fastq",
		host_taxid="2 10239 40674"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml"
	threads: 4
	benchmark:
		dirs_dict["BENCHMARKS"] +"/kraken/{sample}_contaminant_removal.tsv"
	resources:
		mem_gb=40
	shell:
		"""
		python {input.kraken_tools}/extract_kraken_reads.py -k {input.kraken_output_paired} \
			-s1 {input.forward_paired} -s2 {input.reverse_paired} \
			-o {output.forward_paired} -o2 {output.reverse_paired} \
			--exclude --taxid {params.host_taxid} --include-children -r {input.kraken_report_paired} --fastq-output
		python {input.kraken_tools}/extract_kraken_reads.py -k {input.kraken_output_unpaired} \
			-s {input.merged_unpaired} -o {output.unpaired} --exclude --taxid {params.host_taxid} --include-children \
			-r {input.kraken_report_unpaired} --fastq-output
		"""

rule contaminants_KRAKEN_clean:
	input:
		forward_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.fastq.gz"),
		reverse_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.fastq.gz"),
		unpaired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_unpaired_clean.fastq.gz",
		kraken_db=(config['kraken_db']),
		kraken_tools=(config['kraken_tools']),
	output:
		kraken_output_paired=temp(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_output_paired_clean.csv"),
		kraken_report_paired=(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_clean.csv"),
	params:
		kraken_db=config['kraken_db'],
	message:
		"Assesing taxonomy with kraken2 on clean reads"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/kraken/{sample}_clean.tsv"
	priority: 1
	threads: 8
	shell:
		"""
		kraken2 --db {params.kraken_db} --threads {threads} \
			--paired {input.forward_paired} {input.reverse_paired} \
			--output {output.kraken_output_paired} --report {output.kraken_report_paired} \
			--report-minimizer-data
		"""

rule prekrakenMultiQC:
	input:
		expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired.csv", sample=SAMPLES),
	output:
		multiqc=dirs_dict["QC_DIR"]+ "/pre_decontamination_kraken_multiqc_report.html"
		# 		multiqc=dirs_dict["QC_DIR"]+ "/preQC_illumina_report.html",
	params:
		fastqc_dir=dirs_dict["CLEAN_DATA_DIR"],
		html_name="pre_decontamination_kraken_multiqc_report.html",
		multiqc_dir=dirs_dict["QC_DIR"]
	message:
		"Generating MultiQC report kraken (pre QC)"
	priority: 1
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/multiQC/multiqc_kraken_pre.tsv"
	shell:
		"""
		multiqc -f {input} -o {params.multiqc_dir} -n {params.html_name}
		"""

rule postkrakenMultiQC:
	input:
		expand(dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_kraken2_report_paired_clean.csv", sample=SAMPLES),
	output:
		multiqc=dirs_dict["QC_DIR"]+ "/post_decontamination_kraken_multiqc_report.html"
	params:
		fastqc_dir=dirs_dict["CLEAN_DATA_DIR"],
		html_name="post_decontamination_kraken_multiqc_report.html",
		multiqc_dir=dirs_dict["QC_DIR"]
	message:
		"Generating MultiQC report kraken (post QC)"
	priority: 1
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/multiQC/multiqc_kraken_post.tsv"
	shell:
		"""
		multiqc -f {input} -o {params.multiqc_dir} -n {params.html_name}
		"""
		
rule superDeduper_pcr:
	input:
		forward_file=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['forward_tag']) + ".fastq.gz",
		reverse_file=dirs_dict["RAW_DATA_DIR"] + "/{sample}_" + str(config['reverse_tag']) + ".fastq.gz",
	output:
		duplicate_stats=(dirs_dict["QC_DIR"] + "/{sample}_stats_pcr_duplicates.log"),
		deduplicate=temp(dirs_dict["QC_DIR"] + "/{sample}_stats_pcr_duplicates.out"),
	message:
		"Detect PCR duplicates"
	conda:
		dirs_dict["ENVS_DIR"]+ "/QC.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] +"/SuperDeduper/{sample}_pcr_duplicates.tsv"
	shell:
		"""
		hts_SuperDeduper -L {output.duplicate_stats} -1 {input.forward_file} -2 {input.reverse_file} > {output.deduplicate}
		"""
