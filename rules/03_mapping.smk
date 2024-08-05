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
	benchmark:
		dirs_dict["BENCHMARKS"] +"/03_mapping/download_reference_{reference_genome}.tsv"
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml",
	shell:
		"""
		echo {params.reference_genome} > {output.txt_acc}
		bit-dl-ncbi-assemblies -w {output.txt_acc} -f fasta
		gunzip {params.fasta_gz}
		mv {params.fasta_temp} {output.reference_fasta}
		"""

rule genome_Index_fasta:
	input:
		genome="{reference_genome}.fasta",
	output:
		index="{reference_genome}.1.ht2",
	params:
		reference_genome="{reference_genome}" 
	benchmark:
		dirs_dict["BENCHMARKS"] +"/03_mapping/index_reference_{reference_genome}.tsv"
	message:
		"Indexing genome with HISAT2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	resources:
		cpus_per_task= 8,
	shell:
		"""
		hisat2-build -p {resources.cpus_per_task} {input.genome} {params.reference_genome}
		"""

rule map_reads_to_index:
	input:
		index=dirs_dict["GENOMES_DIR"] +"/{reference_genome}.1.ht2", 
		forward_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.fastq.gz",
		reverse_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.fastq.gz",
	output:
		sam=dirs_dict["MAPPING_DIR"] + "/{sample}_{reference_genome}.sam",
		bam=dirs_dict["MAPPING_DIR"] + "/{sample}_{reference_genome}.bam",
		sorted_bam=dirs_dict["MAPPING_DIR"] + "/{sample}_{reference_genome}_sorted.bam",
		flagstats=dirs_dict["MAPPING_DIR"] + "/{sample}_{reference_genome}_flagstats.txt",
	benchmark:
		dirs_dict["BENCHMARKS"] +"/03_mapping/map_reads_to_index_{sample}_{reference_genome}.tsv"
	wildcard_constraints:
		reference_genome="[A-Z]+_[0-9]+.[0-9]"
	params:
		genome=dirs_dict["GENOMES_DIR"] +"/{reference_genome}",
	message:
		"Mapping reads to reference genomes with HISAT2"
	conda:
		dirs_dict["ENVS_DIR"] + "/env1.yaml"
	resources:
		cpus_per_task= 16,
	shell:
		"""
		hisat2 --dta -p {resources.cpus_per_task} -x {params.genome} -1 {input.forward_paired} -2 {input.reverse_paired} -S {output.sam}
      	samtools view  -@ {resources.cpus_per_task} -bS {output.sam}  > {output.bam} 
		samtools sort -@ {resources.cpus_per_task} {output.bam} -o {output.sorted_bam}
		samtools index {output.sorted_bam}
		samtools flagstat {output.sorted_bam} > {output.flagstats}
		"""