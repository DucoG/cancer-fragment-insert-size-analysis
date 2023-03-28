import sys
import pandas as pd
import numpy as np
from scipy import stats


def calculate_summary_statistics(metrics_files, output_file):
    summary_stats = []
    for file in metrics_files:
        data = pd.read_csv(file, sep='\t', names=["insert_size", "count"])
        sample = file.split('/')[-1][:-12]
        mean_insert_size = np.average(data["insert_size"], weights=data["count"])
        median_insert_size = np.average(data["insert_size"], weights=data["count"])

        # calculate mode by finding the largest coubnt and then finding the insert size with that count
        mode_insert_size = data.loc[data["count"] == data["count"].max(), "insert_size"].iloc[0]

        std_dev_insert_size = np.sqrt(np.average((data["insert_size"]-mean_insert_size)**2, weights=data["count"]))
        summary_stats.append([sample, mean_insert_size, median_insert_size, mode_insert_size, std_dev_insert_size])

    summary_df = pd.DataFrame(summary_stats, columns=["sample", "mean", "median", "mode", "std_dev"])
    summary_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    metrics_files = sys.argv[1:-1]
    output_file = sys.argv[-1]
    calculate_summary_statistics(metrics_files, output_file)