rule download_reference_genomes:
	output:
		reference_fasta=dirs_dict["GENOMES_DIR"] +"/{reference_genome}.fasta",
		txt_acc=temp((dirs_dict["GENOMES_DIR"] +"/{reference_genome}.txt")),
	message:
		"Downloading reference genomes"
	params:
		reference_genome="{reference_genome}",
		fasta_gz=("{reference_genome}.fa.gz"),
		fasta_temp=("{reference_genome}.fa"),
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml",
	threads: 1
	shell:
		"""
		echo {params.reference_genome} > {output.txt_acc}
		bit-dl-ncbi-assemblies -w {output.txt_acc} -f fasta
		gunzip {output.fasta_temp}
		mv {params.fasta_temp} {output.reference_fasta}
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
		sam=dirs_dict["MAPPING_DIR"] + "/{sample}_{genome_name}.sam",
	wildcard_constraints:
		genome_name="[A-Z]+_[0-9]+.[0-9]"
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