import sys
import pandas as pd
import matplotlib.pyplot as plt

def plot_histogram(metrics_file, output_file):
    data = pd.read_csv(metrics_file, sep='\t', names=["insert_size", "count"])
    plt.figure(figsize=(10, 6))
    plt.bar(data["insert_size"], data["count"], width=1)
    plt.xlabel("Insert Size")
    plt.ylabel("Count")
    plt.title("Fragment Length Distribution")
    plt.savefig(output_file, dpi=300)
    plt.close()

if __name__ == "__main__":
    metrics_file = sys.argv[1]
    output_file = sys.argv[2]
    plot_histogram(metrics_file, output_file)
