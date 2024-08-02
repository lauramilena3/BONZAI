rule download_reference_genomes:
	output:
		contaminant_fasta=dirs_dict["GENOMES_DIR"] +"/{reference_genome}.fasta",
		contaminant_dir=temp(directory(dirs_dict["GENOMES_DIR"] +"/temp_{reference_genome}")),
	message:
		"Downloading reference genomes"
	params:
		contaminants_dir=dirs_dict["GENOMES_DIR"],
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml",
	threads: 16
	shell:
		"""
		mkdir {output.contaminant_dir}
		cd {output.contaminant_dir}
		wget $(esearch -db "assembly" -query {wildcards.contaminant} | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F"/" '{{print $0"/"$NF"_genomic.fna.gz"}}')
		gunzip -f *gz
		cat *fna >> {output.contaminant_fasta}
		"""

rule genome_Index_fasta:
	input:
		genome="{genome_name}.fasta",
	output:
		index="{genome_name}.1.ht2",
	params:
		genome_name="{genome_name}" 
	message:
		"Indexing genome with HISAT2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 8
	shell:
		"""
		hisat2-build -p {threads} {input.genome} {params.genome_name}
		"""

rule map_reads_to_index:
	input:
		index=dirs_dict["GENOMES_DIR"] +"/{genome_name}.1.ht2", 
		forward_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.fastq.gz",
		reverse_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.fastq.gz",
	output:
		sam=dirs_dict["MAPPING_DIR"] + "/{genome_name}_{sample}.sam",
	params:
		genome="{genome_name}" 
	message:
		"Indexing genome with HISAT2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	threads: 16
	shell:
		"""
		hisat2 --dta -p {threads} -x {params.genome} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam}
		"""