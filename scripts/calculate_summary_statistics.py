import sys
import pandas as pd

def calculate_summary_statistics(metrics_files, output_file):
    summary_stats = []
    for file in metrics_files:
        data = pd.read_csv(file, sep='\t', names=["insert_size", "count"])
        sample = file.split('/')[-1][:-12]
        mean = data["count"].mean()
        median = data["count"].median()
        mode = data["count"].mode().get(0, None)
        std_dev = data["count"].std()
        summary_stats.append([sample, mean, median, mode, std_dev])

    summary_df = pd.DataFrame(summary_stats, columns=["sample", "mean", "median", "mode", "std_dev"])
    summary_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    metrics_files = sys.argv[1:-1]
    output_file = sys.argv[-1]
    calculate_summary_statistics(metrics_files, output_file)