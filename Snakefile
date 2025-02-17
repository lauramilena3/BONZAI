import os
import re
import glob
import pandas as pd

#======================================================
# Config files
#======================================================
configfile: "config.yaml"

#======================================================
# Global variables
#======================================================

RAW_DATA_DIR = config['input_dir'].rstrip("/")
RESULTS_DIR = config['results_dir'].rstrip("/")

if not RESULTS_DIR and RAW_DATA_DIR:
    RESULTS_DIR = os.path.abspath(os.path.join(RAW_DATA_DIR, os.pardir))

PAIRED = True

if RAW_DATA_DIR:
    sample_files = glob.glob(RAW_DATA_DIR + '/*_L00*_R[12]_001.fastq.gz')
    SAMPLES = sorted(set(re.sub(r'_L00[0-9]_R[12]_001.fastq.gz$', '', os.path.basename(f)) for f in sample_files))
else:
    RAW_DATA_DIR = RESULTS_DIR + "/00_RAW_DATA"
    SAMPLES = []

# Verificar si hay muestras detectadas
if not SAMPLES:
    print("⚠️ No samples found! Please check your input directory.")

# Definir directorios

dir_list = [
    "RULES_DIR", "ENVS_DIR", "DB_DIR", "ADAPTERS_DIR", "RAW_NOTEBOOKS", "RAW_DATA_DIR",
    "QC_DIR", "CLEAN_DATA_DIR", "ASSEMBLY_DIR", "MAPPING_DIR", "TRANSCRIPTOME_ASSEMBLY_DIR",
    "ANNOTATION", "ABUNDANCE", "BENCHMARKS", "NOTEBOOKS_DIR", "PLOTS_DIR", "GENOMES_DIR",
    "TEMP_CLUSTER_DIR"
]

dir_names = [
    "rules", "../envs", "db", "db/adapters", "../notebooks", RAW_DATA_DIR,
    RESULTS_DIR + "/01_QC", RESULTS_DIR + "/01_QC", RESULTS_DIR + "/02_ASSEMBLY",
    RESULTS_DIR + "/03_MAPPING", RESULTS_DIR + "/04_TRANSCRIPTOME_ASSEMBLY",
    RESULTS_DIR + "/05_ANNOTATION", RESULTS_DIR + "/06_GENE_EXPRESSION",
    RESULTS_DIR + "/BENCHMARK", RESULTS_DIR + "/NOTEBOOKS", RESULTS_DIR + "/FIGURES_AND_TABLES",
    RESULTS_DIR + "/REFERENCE_GENOMES", RESULTS_DIR + "/TEMP_CLUSTER_DIR"
]

dirs_dict = dict(zip(dir_list, dir_names))

READ_TYPES = [config['forward_tag'], config['reverse_tag']]
REFERENCE_GENOME_ACC = config['reference_acc_list'].split()

print("📂 Input Directory:", RAW_DATA_DIR)
print("📂 Results Directory:", RESULTS_DIR)
print("📂 Sample Names:", ", ".join(SAMPLES) if SAMPLES else "No samples found")
print("📂 Read Types:", READ_TYPES)
print("📂 Reference Genomes:", REFERENCE_GENOME_ACC)

# ----------------------------------------------------------------------------
# CREAR LOS DIRECTORIOS SI NO EXISTEN
# ----------------------------------------------------------------------------

for dir_path in dirs_dict.values():
    os.makedirs(dir_path, exist_ok=True)

print("\n🔹 Se han creado/verificado los siguientes directorios:")
for key, value in dirs_dict.items():
    print(f"🔸 {key}: {value}")

# ----------------------------------------------------------------------------
# FUNCIONES PARA OBTENER LOS ARCHIVOS DE ENTRADA
# ----------------------------------------------------------------------------

def inputReadsCount(wildcards):
    if not SAMPLES:
        return []
    return expand(dirs_dict["RAW_DATA_DIR"] + "/{sample}_L00{lane}_R{read}_read_count.txt", 
                  sample=SAMPLES, lane=["1", "2", "3", "4"], read=READ_TYPES)

def inputQC(wildcards):
    return [
        dirs_dict["QC_DIR"] + "/preQC_illumina_report.html",
        dirs_dict["QC_DIR"] + "/postQC_illumina_report.html",
        dirs_dict["QC_DIR"] + "/pre_decontamination_kraken_multiqc_report.html",
        dirs_dict["QC_DIR"] + "/post_decontamination_kraken_multiqc_report.html"
    ]

def inputMapping(wildcards):
    if not SAMPLES or not REFERENCE_GENOME_ACC:
        return []
    return expand(dirs_dict["MAPPING_DIR"] + "/{sample}_{reference_genome}.sam", 
                  sample=SAMPLES, reference_genome=REFERENCE_GENOME_ACC)

def inputDeNovoAssembly(wildcards):
    if not SAMPLES:
        return []
    return expand(dirs_dict["ASSEMBLY_DIR"] + "/{sample}_trinity/{sample}.fasta", 
                  sample=SAMPLES)

# ----------------------------------------------------------------------------
# REGLAS DEL PIPELINE
# ----------------------------------------------------------------------------

rule all:
    input:
        inputReadsCount({}) if SAMPLES else [],
        inputQC({}),
        # inputDeNovoAssembly({}) if SAMPLES else [],
        # inputMapping({}) if SAMPLES and REFERENCE_GENOME_ACC else [],
        # inputTranscriptomeAssembly({}),
        # inputAnnotation({}),
        # inputAbundance({}),

include: os.path.join(dirs_dict["RULES_DIR"], '00_download_tools.smk')
include: os.path.join(dirs_dict["RULES_DIR"], '01_quality_control.smk')
include: os.path.join(dirs_dict["RULES_DIR"], '02_de_novo_assembly.smk')
include: os.path.join(dirs_dict["RULES_DIR"], '03_mapping.smk')
include: os.path.join(dirs_dict["RULES_DIR"], '04_transcriptome_assembly.smk')
# include: os.path.join(dirs_dict["RULES_DIR"], '05_annotation.smk')
# include: os.path.join(dirs_dict["RULES_DIR"], '06_gene_expression.smk')