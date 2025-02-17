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

RAW_DATA_DIR = config["dirs"]["RAW_DATA_DIR"].rstrip("/")
RESULTS_DIR = os.path.dirname(RAW_DATA_DIR)

PAIRED = True

# Asegurar que los READ_TYPES sean correctos
READ_TYPES = [config["forward_tag"], config["reverse_tag"]]
if any(rt not in ["R1", "R2"] for rt in READ_TYPES):
    raise ValueError(f"⚠️ Error en READ_TYPES: {READ_TYPES}. Debe ser ['R1', 'R2'].")

# Buscar muestras en el directorio de datos crudos
if RAW_DATA_DIR:
    sample_files = glob.glob(RAW_DATA_DIR + '/*_L00*_R[12]_001.fastq.gz')
    SAMPLES = sorted(set(re.sub(r'_L00\d+_R[12]_001.fastq.gz$', '', os.path.basename(f)) for f in sample_files))
else:
    RAW_DATA_DIR = RESULTS_DIR + "/00_RAW_DATA"
    SAMPLES = []

# Verificar si hay muestras detectadas
if not SAMPLES:
    print("⚠️ No samples found! Please check your input directory.")

REFERENCE_GENOME_ACC = config["reference_acc_list"].split()

print("📂 Input Directory:", RAW_DATA_DIR)
print("📂 Results Directory:", RESULTS_DIR)
print("📂 Sample Names:", ", ".join(SAMPLES) if SAMPLES else "No samples found")
print("📂 Read Types:", READ_TYPES)
print("📂 Reference Genomes:", REFERENCE_GENOME_ACC)

# ----------------------------------------------------------------------------
# CREAR LOS DIRECTORIOS SI NO EXISTEN
# ----------------------------------------------------------------------------

for key, dir_path in config["dirs"].items():
    os.makedirs(dir_path, exist_ok=True)

print("\n🔹 Se han creado/verificado los siguientes directorios:")
for key, value in config["dirs"].items():
    print(f"🔸 {key}: {value}")

# ----------------------------------------------------------------------------
# FUNCIONES PARA OBTENER LOS ARCHIVOS DE ENTRADA
# ----------------------------------------------------------------------------

def inputReadsCount(wildcards):
    if not SAMPLES:
        return []
    return expand(config["dirs"]["RAW_DATA_DIR"] + "/{sample}_L00{lane}_R{read}_read_count.txt", 
                  sample=SAMPLES, lane=["1", "2", "3", "4"], read=READ_TYPES)

def inputQC(wildcards):
    return [
        config["dirs"]["QC_DIR"] + "/preQC_illumina_report.html",
        config["dirs"]["QC_DIR"] + "/postQC_illumina_report.html",
        config["dirs"]["QC_DIR"] + "/pre_decontamination_kraken_multiqc_report.html",
        config["dirs"]["QC_DIR"] + "/post_decontamination_kraken_multiqc_report.html"
    ]

def inputMapping(wildcards):
    if not SAMPLES or not REFERENCE_GENOME_ACC:
        return []
    return expand(config["dirs"]["MAPPING_DIR"] + "/{sample}_{reference_genome}.sam", 
                  sample=SAMPLES, reference_genome=REFERENCE_GENOME_ACC)

def inputDeNovoAssembly(wildcards):
    if not SAMPLES:
        return []
    return expand(config["dirs"]["ASSEMBLY_DIR"] + "/{sample}_trinity/{sample}.fasta", 
                  sample=SAMPLES)

# ----------------------------------------------------------------------------
# REGLAS DEL PIPELINE
# ----------------------------------------------------------------------------

rule all:
    input:
        inputReadsCount({}) if SAMPLES else [],
        inputQC({}),

include: os.path.join("rules", '00_download_tools.smk')
include: os.path.join("rules", '01_quality_control.smk')
include: os.path.join("rules", '02_de_novo_assembly.smk')
include: os.path.join("rules", '03_mapping.smk')
include: os.path.join("rules", '04_transcriptome_assembly.smk')