#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot the chocolate side for each algorithm.")
    parser.add_argument('--dataset', required=True, help='Dataset name.')
    parser.add_argument('--summary_stats', nargs='+', required=True, help='Summary statistics CSV files.')
    parser.add_argument('--qps_files', nargs='+', required=True, help='QPS CSV files.')
    parser.add_argument('--output', required=True, help='Output plot file.')
    return parser.parse_args()

def find_best_parameters(summary_stats_files):
    best_params = {}
    for file in summary_stats_files:
        # Extract algorithm name from filename
        basename = os.path.basename(file)
        parts = basename.split('_')
        dataset = parts[0]
        algorithm = parts[1]

        df = pd.read_csv(file)
        if 'AUC' not in df.columns:
            continue

        # Find the row with the highest AUC
        max_auc = df['AUC'].max()
        best_row = df[df['AUC'] == max_auc].iloc[0]

        # Store the parameters
        best_params[algorithm] = {
            'params': best_row.to_dict(),
            'summary_file': file
        }
    return best_params

def get_qps_file_for_best_param(qps_files, algorithm, params):
    # Assuming QPS files have parameters encoded in their filenames similar to summary_stats files
    for file in qps_files:
        if algorithm in file:
            # Simple matching based on filename; adjust as needed
            if all(f"_{key}{value}" in file for key, value in params.items() if key != 'AUC'):
                return file
    return None

def plot_chocolate_side(best_params, qps_files, dataset, output_path):
    plt.figure(figsize=(10, 6))

    for algorithm, data in best_params.items():
        params = data['params']
        qps_file = get_qps_file_for_best_param(qps_files, algorithm, params)

        if not qps_file:
            print(f"No QPS file found for algorithm {algorithm} with parameters {params}")
            continue

        qps_df = pd.read_csv(qps_file)
        if 'recall' not in qps_df.columns or 'QPS' not in qps_df.columns:
            print(f"QPS file {qps_file} does not contain 'recall' or 'QPS' columns.")
            continue

        plt.plot(qps_df['recall'], qps_df['QPS'], marker='o', label=f"{algorithm}")

    plt.xlabel('Recall')
    plt.ylabel('QPS')
    plt.title(f'Chocolate Side Plot for Dataset: {dataset}')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Chocolate side plot saved to {output_path}")

def main():
    args = parse_arguments()

    # Find the best parameters for each algorithm
    best_params = find_best_parameters(args.summary_stats)

    # Plot the chocolate side
    plot_chocolate_side(best_params, args.qps_files, args.dataset, args.output)

if __name__ == "__main__":
    main()
