# Ovarian Cancer Paired-end Sequencing Insert Size Analysis

This repository contains a Snakemake pipeline and accompanying scripts for analyzing shallow whole-genome sequencing (sWGS) data from ctDNA samples to differentiate between benign and malignant ovarian cancer cases using fragment length distribution.

## Overview

The pipeline performs the following steps:

1. Calculate insert size metrics from BAM files for each sample.
2. Generate fragment length distribution histograms for each sample.
3. Generate a combined histogram comparing mean insert size distributions, including standard deviation, for benign and malignant cases.
4. Calculate summary statistics (mean, median, mode, and standard deviation) for each case.

## Requirements

- Python 3.9
- Snakemake
- Samtools
- Conda (optional, for managing environments)

## Usage

1. Clone this repository to your local machine:

```
git clone https://github.com/DucoG/cancer_insert_size_analysis.git
```

2. Modify the `config_file.yaml` to include the paths to your data files (BAM and BAI files) and the Conda environment file, if using Conda.

3. Install the Conda environment:

```
conda env create -f env.yaml
```

4. Activate the Conda environment:

```
conda activate wgs
```

5. Run the Snakemake pipeline with the `--use-conda` argument:

```
snakemake --cores N --use-conda
```

Replace `N` with the number of cores you want to use for the analysis.

## Output

The pipeline generates the following output files:

- Insert size metrics for each sample: `insert_size_metrics/`
- Fragment length distribution histograms for each sample: `insert_size_histograms/`
- A combined histogram comparing mean insert size distributions for benign and malignant cases: `combined_histogram/combined_histogram.pdf`
- Plotting data (mean and standard deviation) for the combined histogram: `combined_histogram/plotting_data.csv`
- Summary statistics (mean, median, mode, and standard deviation) for each case: `summary_statistics/summary_statistics.csv`

## License

[MIT License](LICENSE)