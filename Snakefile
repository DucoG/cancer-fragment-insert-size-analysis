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
        expand("nucleosome_distance/{sample}_nucleosome_distance.json", sample=SAMPLES)


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

# rule to generate a json file for each sample containing the read count per distance to a nucleosome and the flagcount for each read. 
# (wgs) d.gaillard@darwin:~/paired_ovarian/fragment_lengh_distibution$ python scripts/paired_calculate_nucleosome_distance.py --help
# Usage: paired_calculate_nucleosome_distance.py [OPTIONS] BAMFILE_PATH
#                                                OUTPUT_PATH

#   Extracts distance to closest nucleosome for each read from a bamfile.

#   Args:

#       bamfile_path (str): path of input file. If it doesnt end with .bam, it
#       is assumed to be a text file containing the path for a bam file on each
#       line.

#       output_path (str): path to output the resulting count array

#       nulceosome_list (str): path to list of nucleosome locations per contig
#       in tsv format

#       unwanted_chrs (str): unwanted chr

# Options:
#   --nucleosome_list PATH    path to list containing nucleosome locations on
#                             contigs in a tsv format  [required]
#   -u, --unwanted_chrs TEXT  chromosomes to leave out in the analysis.
#                             Structure to be used: chr1
#   --help                    Show this message and exit.

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
        "python scripts/paired_calculate_nucleosome_distance.py {input.bam} {output.json} --nucleosome_list {input.nuc_list} -u chrY"

