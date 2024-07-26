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
		kraken_db=directory(config['kraken_db']),
	message:
		"Downloading miniKraken database"
	threads: 4
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	shell:
		"""
		wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz
		mkdir {output.kraken_db}
		tar -
		"""
