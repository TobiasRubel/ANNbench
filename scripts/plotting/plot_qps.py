import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import warnings

def plot_qps_vs_recall(csv_files, plot_file):
    plt.figure(figsize=(10, 6))

    for csv_file in csv_files:
        data = pd.read_csv(csv_file)
        # check if data is "No data extracted from log file"
        if data.shape[1] <= 1:
            warnings.warn(f"No data extracted from {csv_file}")
            continue
        algorithm = os.path.basename(csv_file).split('_')[1:-2]  # Assuming algorithm is in the file name
        
        plt.plot(data['recall'], data['QPS'], marker='o', label=algorithm)

    k = os.path.basename(csv_files[0]).split('_')[-2][1:]
    nthreads = os.path.basename(csv_files[0]).split('_')[-1].split('.')[0][8:]
    plt.xlabel('Recall')
    plt.ylabel('QPS')
    plt.title('QPS vs. Recall ({}@{}, nthreads = {})'.format(k,k,nthreads))
    plt.legend()

    plt.savefig(plot_file)
    plt.close()
    print(f"Plot saved to {plot_file}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python plot_qps_vs_recall.py <csv_files> <plot_file>")
        sys.exit(1)

    csv_files = sys.argv[1:-1]
    plot_file = sys.argv[-1]

    plot_qps_vs_recall(csv_files, plot_file)
