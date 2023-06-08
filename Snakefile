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
        "summary_statistics/fragment_length_statistics.csv",
        "summary_statistics/mito_ratios.csv",
        expand("nucleosome_distance/{sample}_nucleosome_distance.json", sample=SAMPLES),
        expand('binwise_fragmentomics/{sample}_b1000000_fragmentomics.csv', sample=SAMPLES),
        expand('binwise_fragmentomics/{sample}_b5000000_fragmentomics.csv', sample=SAMPLES),
        expand('binwise_fragmentomics/{sample}_b500000_fragmentomics.csv', sample=SAMPLES),
        expand('binwise_fragmentomics/{sample}_b100000_fragmentomics.csv', sample=SAMPLES),
        expand('binwise_fragmentomics/{sample}_b50000_fragmentomics.csv', sample=SAMPLES),
        expand('binwise_fragmentomics/{sample}_b250000_fragmentomics.csv', sample=SAMPLES),
        expand('delfi_ratios/{sample}_b{binsize}_delfi_ratios.csv', sample=SAMPLES, binsize=[5000000, 1000000, 500000, 250000])
        





rule samtools_stats:
    input:
        bam=lambda wildcards: BAM_PATHS[wildcards.sample],
        bai=lambda wildcards: BAM_PATHS[wildcards.sample] + ".bai"
    output:
        stats="samtools_stats/{sample}_stats.txt"
    conda:
        config['conda_env']
    shell:
        "samtools stats {input.bam} > {output.stats}"

rule calculate_insert_size_metrics:
    input:
        stats="samtools_stats/{sample}_stats.txt"
    output:
        metrics="insert_size_metrics/{sample}_metrics.txt"
    conda:
        config['conda_env']
    shell:
        "grep '^IS' {input.stats} | "
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
        summary_csv="summary_statistics/fragment_length_statistics.csv"
    conda:
        config['conda_env']
    shell:
        "python scripts/calculate_summary_statistics.py {input} {output.summary_csv}"

# rule to generate a csv file containing the ratio of mitochondrial reads to total reads, takes in all bam files and outputs a csv file with the sample name and the ratio in the columns
rule calculate_mito_ratios:
    input:
    # use expand to generate a list of all the bam files
        sample_info["path"]
    output:
        "summary_statistics/mito_ratios.csv"
    conda:
        config['conda_env']
    shell:
        "python scripts/calculate_mito_ratios.py {output} {input}"


rule paired_calculate_nucleosome_distance:
    input:
        bam=lambda wildcards: BAM_PATHS[wildcards.sample],
        bai=lambda wildcards: BAM_PATHS[wildcards.sample] + ".bai",
        nuc_list="data/nuc_center_list.txt"
    output:
        json="nucleosome_distance/{sample}_nucleosome_distance.json"
    conda:
        config['conda_env']
    shell:
        "python scripts/paired_calculate_nucleosome_distance.py {input.bam} {output.json} --nucleosome_list {input.nuc_list} -u chrY -u chrM"

rule calculate_binwise_fragmentomics:
    input:
        bam=lambda wildcards: BAM_PATHS[wildcards.sample],
        bai=lambda wildcards: BAM_PATHS[wildcards.sample] + ".bai",
        nuc_list="data/nuc_center_list.txt"
    output:
        "binwise_fragmentomics/{sample}_b{binsize}_fragmentomics.csv"
    conda:
        config['conda_env']
    shell:
        "python scripts/binwise_fragmentomics_analyzer.py {input.bam} {input.nuc_list} {wildcards.binsize} chrM {output}"

rule calculate_delfi_ratios:
    input:
        bam=lambda wildcards: BAM_PATHS[wildcards.sample],
        bai=lambda wildcards: BAM_PATHS[wildcards.sample] + ".bai",
        fasta="/home/d.gaillard/source/PEsWGS-alignment-snakemake/ref_genome/hg19.fa"
    output:
        delfi_out="delfi_ratios/{sample}_b{binsize}_delfi_ratios.csv",
        fragment_distributions="binwise_fragment_length_distributions/{sample}_b{binsize}_fragment_length_distribution.json"
    conda:
        config['conda_env']
    shell:
        "python scripts/calculate_delfi_ratios.py {input.bam} {input.fasta} {wildcards.binsize} {output.delfi_out} {output.fragment_distributions} chrM"