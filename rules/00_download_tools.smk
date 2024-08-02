rule getKrakenTools:
	output:
		kraken_tools=directory(config['kraken_tools']),
	message:
		"Downloading KrakenTools"
	threads: 4
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	shell:
		"""
		mkdir -p tools
		cd tools
		git clone https://github.com/jenniferlu717/KrakenTools
		chmod 777 KrakenTools/*
		"""

rule downloadminiKrakenDB:
	output:
		minikraken_db=directory("db/KRAKEN/minikraken_8GB_20200312"),
		minikraken_gz=temp("minikraken_8GB_202003.tgz"),
	params:
		db_dir="db/KRAKEN/"
	message:
		"Downloading miniKraken database"
	threads: 4
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	shell:
		"""
		wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz
		mkdir -p {params.db_dir}
		tar -xvf {output.minikraken_gz} -C {params.db_dir}
		"""

rule download_reference_genomes:
	output:
		contaminant_fasta=dirs_dict["GENOMES_DIR"] +"/{contaminant}.fasta",
		contaminant_dir=temp(directory(dirs_dict["GENOMES_DIR"] +"/temp_{contaminant}")),
	message:
		"Downloading reference genomes"
	params:
		contaminants_dir=dirs_dict["GENOMES_DIR"],
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1.yaml",
	threads:
		16
	shell:
		"""
		mkdir {output.contaminant_dir}
		cd {output.contaminant_dir}
		wget $(esearch -db "assembly" -query {wildcards.contaminant} | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F"/" '{{print $0"/"$NF"_genomic.fna.gz"}}')
		gunzip -f *gz
		cat *fna >> {output.contaminant_fasta}
		"""