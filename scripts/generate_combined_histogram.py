import sys
import pandas as pd
import matplotlib.pyplot as plt

def calculate_mean_distribution(metrics_files):
    combined_data = pd.concat([pd.read_csv(file, sep='\t', names=["insert_size", "count"]) for file in metrics_files])
    mean_data = combined_data.groupby("insert_size").mean().reset_index()
    std_data = combined_data.groupby("insert_size").std().reset_index()
    return mean_data, std_data

def plot_combined_histogram(benign_metrics_files, malignant_metrics_files, output_file, plotting_data_file):
    benign_mean, benign_std = calculate_mean_distribution(benign_metrics_files)
    malignant_mean, malignant_std = calculate_mean_distribution(malignant_metrics_files)

    plt.figure(figsize=(10, 6))

    plt.plot(benign_mean["insert_size"], benign_mean["count"], color='blue', label='Benign')
    plt.fill_between(benign_mean["insert_size"], benign_mean["count"] - benign_std["count"], benign_mean["count"] + benign_std["count"], color='blue', alpha=0.1)

    plt.plot(malignant_mean["insert_size"], malignant_mean["count"], color='red', label='Malignant')
    plt.fill_between(malignant_mean["insert_size"], malignant_mean["count"] - malignant_std["count"], malignant_mean["count"] + malignant_std["count"], color='red', alpha=0.1)

    plt.xlabel("Insert Size")
    plt.ylabel("Mean Count")
    plt.title("Fragment Length Distribution - Benign vs Malignant")
    plt.legend()
    plt.savefig(output_file, dpi=300)
    plt.close()

    # Save plotting data to a file
    benign_mean.columns = ["insert_size", "benign_mean"]
    benign_std.columns = ["insert_size", "benign_std"]
    malignant_mean.columns = ["insert_size", "malignant_mean"]
    malignant_std.columns = ["insert_size", "malignant_std"]
    
    plotting_data = pd.merge(benign_mean, benign_std, on="insert_size")
    plotting_data = pd.merge(plotting_data, malignant_mean, on="insert_size")
    plotting_data = pd.merge(plotting_data, malignant_std, on="insert_size")
    
    plotting_data.to_csv(plotting_data_file, index=False)

if __name__ == "__main__":
    try:
        double_hyphen_index = sys.argv.index("--")
    except ValueError:
        raise ValueError("Please separate benign and malignant metrics files with a '--' in the command line arguments.")
    
    benign_metrics_files = sys.argv[1:double_hyphen_index]
    malignant_metrics_files = sys.argv[double_hyphen_index+1:-2]
    output_file = sys.argv[-2]
    plotting_data_file = sys.argv[-1]
    plot_combined_histogram(benign_metrics_files, malignant_metrics_files, output_file, plotting_data_file)