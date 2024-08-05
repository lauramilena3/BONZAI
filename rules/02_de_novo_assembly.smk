rule de_novo_assembly:
	input:
		forward_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_forward_paired_clean.fastq.gz",
		reverse_paired=dirs_dict["CLEAN_DATA_DIR"] + "/{sample}_reverse_paired_clean.fastq.gz",
	output:
		trinity_fasta_temp=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_trinity/Trinity.fasta",
		trinity_fasta=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_trinity/{sample}.fasta",
		jellyfish_count=temp(dirs_dict["TEMP_CLUSTER_DIR"] +"/{sample}_mer_counts.jf")
	params:
		trinity_dir=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_trinity",
	conda:
		dirs_dict["ENVS_DIR"] + "/env2.yaml"
	threads: 16
	benchmark:
		dirs_dict["BENCHMARKS"] +"/deNovoAssembly/{sample}_denovo_assembly.tsv"
	resources:
		mem_gb=100,
		runtime= 60, # minutes
	shell:
		"""
		# jellyfish count -o {output.jellyfish_count}
		Trinity --seqType fq --left {input.forward_paired}  --right {input.reverse_paired}  --output {params.trinity_dir} --max_memory {resources.mem_gb}G --full_cleanup --CPU {threads}
		"""