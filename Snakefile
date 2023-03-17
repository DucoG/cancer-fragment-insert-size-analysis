import os
import pandas as pd

configfile: "config_file.yaml"

# Read CSV file with paths and labels
sample_info = pd.read_csv(config['datafile'])

# Extract sample names, BAM paths, and labels from the CSV
SAMPLES = [os.path.splitext(os.path.basename(path))[0] for path in sample_info["path"]]
BAM_PATHS = dict(zip(SAMPLES, sample_info["path"]))
LABELS = dict(zip(SAMPLES, sample_info["label"]))

BENIGN_SAMPLES = [s for s in SAMPLES if LABELS[s] == "benign"]
MALIGNANT_SAMPLES = [s for s in SAMPLES if LABELS[s] == "malignant"]

rule all:
    input:
        expand("insert_size_metrics/{sample}_metrics.txt", sample=SAMPLES),
        expand("insert_size_histograms/{sample}_histogram.pdf", sample=SAMPLES),
        "combined_histogram/combined_histogram.pdf",
        "summary_statistics/summary_statistics.csv"


rule calculate_insert_size_metrics:
    input:
        bam=lambda wildcards: BAM_PATHS[wildcards.sample],
        bai=lambda wildcards: BAM_PATHS[wildcards.sample] + ".bai"
    output:
        metrics="insert_size_metrics/{sample}_metrics.txt"
    conda:
        config['conda_env']
    shell:
        "samtools stats {input.bam} | grep '^IS' | "
        "awk '{{print $2 \"\t\" $3}}' > {output.metrics}"

rule generate_histogram:
    input:
        metrics="insert_size_metrics/{sample}_metrics.txt"
    output:
        histogram="insert_size_histograms/{sample}_histogram.pdf"
    conda:
        config['conda_env']
    shell:
        "python scripts/plot_histogram.py {input.metrics} {output.histogram}"

rule generate_combined_histogram:
    input:
        benign_metrics=expand("insert_size_metrics/{sample}_metrics.txt", sample=BENIGN_SAMPLES),
        malignant_metrics=expand("insert_size_metrics/{sample}_metrics.txt", sample=MALIGNANT_SAMPLES)
    output:
        histogram="combined_histogram/combined_histogram.pdf",
        plotting_data="combined_histogram/plotting_data.csv"
    conda:
        config['conda_env']
    shell:
        "python scripts/generate_combined_histogram.py "
        "{input.benign_metrics} -- {input.malignant_metrics} "
        "{output.histogram} {output.plotting_data}"

rule calculate_summary_statistics:
    input:
        expand("insert_size_metrics/{sample}_metrics.txt", sample=SAMPLES)
    output:
        summary_csv="summary_statistics/summary_statistics.csv"
    conda:
        config['conda_env']
    shell:
        "python scripts/calculate_summary_statistics.py {input} {output.summary_csv}"
