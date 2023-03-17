import sys
import pandas as pd
import matplotlib.pyplot as plt

def calculate_mean_distribution(metrics_files):
    # Read each file and store in a list of dataframes
    dataframes = [pd.read_csv(file, sep='\t', names=["insert_size", f"count_{i}"]) for i, file in enumerate(metrics_files)]

    # Merge the dataframes on the insert_size column
    combined_data = dataframes[0]
    for df in dataframes[1:]:
        combined_data = combined_data.merge(df, on="insert_size", how="outer")

    # Fill NaN values with 0
    combined_data.fillna(0, inplace=True)

    # Divide each column by the column sum
    for column in combined_data.columns[1:]:
        combined_data[column] /= combined_data[column].sum()

    # Calculate mean and standard deviation per insert size
    mean_data = combined_data.iloc[:, 1:].mean(axis=1).to_frame(name="mean")

    std_data = combined_data.iloc[:, 1:].std(axis=1).to_frame(name="std")
    
    # Add the insert_size column to the mean and std dataframes
    mean_data.insert(0, "insert_size", combined_data["insert_size"])
    std_data.insert(0, "insert_size", combined_data["insert_size"])

    return mean_data, std_data



def plot_combined_histogram(benign_metrics_files, malignant_metrics_files, output_file, plotting_data_file):
    benign_mean, benign_std = calculate_mean_distribution(benign_metrics_files)
    malignant_mean, malignant_std = calculate_mean_distribution(malignant_metrics_files)

    plt.figure(figsize=(10, 6))

    plt.plot(benign_mean["insert_size"], benign_mean["mean"], color='blue', label='Benign')
    plt.fill_between(benign_mean["insert_size"], benign_mean["mean"] - benign_std["std"], benign_mean["mean"] + benign_std["std"], color='blue', alpha=0.1)

    plt.plot(malignant_mean["insert_size"], malignant_mean["mean"], color='red', label='Malignant')
    plt.fill_between(malignant_mean["insert_size"], malignant_mean["mean"] - malignant_std["std"], malignant_mean["mean"] + malignant_std["std"], color='red', alpha=0.1)

    # add limits
    # plt.xlim([0,400])

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