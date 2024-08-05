rule getKrakenTools:
	output:
		kraken_tools=directory(config['kraken_tools']),
	message:
		"Downloading KrakenTools"
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
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	shell:
		"""
		wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz
		mkdir -p {params.db_dir}
		tar -xvf {output.minikraken_gz} -C {params.db_dir}
		"""

