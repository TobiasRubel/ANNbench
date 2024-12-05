#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot the chocolate side for each algorithm.")
    parser.add_argument('--dataset', required=True, help='Dataset name.')
    parser.add_argument('--summary_stats', nargs='+', required=True, help='Summary statistics CSV files.')
    parser.add_argument('--qps_files', nargs='+', required=True, help='QPS CSV files.')
    parser.add_argument('--output', required=True, help='Output plot file.') 
    return parser.parse_args()

def find_best_parameters(summary_stats_files):
    AUC = 'AUC_qps_vs_recall'
    best_params = {}
    for file in summary_stats_files:
        # Extract dataset name from the file path
        # Assuming the file path is results/{dataset}/summary_stats/{filename}
        dataset = os.path.basename(os.path.dirname(os.path.dirname(file)))
        
        # Extract algorithm name from the filename
        basename = os.path.basename(file)
        
        # Remove the dataset name and the underscore from the beginning of the filename
        if basename.startswith(dataset + '_'):
            remainder = basename[len(dataset) + 1:]  # +1 accounts for the underscore
        else:
            print(f"Warning: Filename {basename} does not start with dataset name {dataset}")
            remainder = basename
        
        # The algorithm name is the next part before the next underscore
        parts = remainder.split('_')
        if len(parts) > 0:
            algorithm = parts[0]
        else:
            print(f"Warning: Unable to extract algorithm name from filename {basename}")
            continue
        
        # Proceed with reading the CSV and finding the best parameters
        df = pd.read_csv(file)
        if AUC not in df.columns:
            print(f"Warning: 'AUC' column not found in file {file}")
            continue
        
        # Find the AUC value
        max_auc = df[AUC].max()

        #get the parameters for this run from the filename
        params = parts[0:-1]

        #update best_params for this algorithm if it has a higher AUC
        try:
            if best_params[algorithm]['AUC'] < max_auc:
                best_params[algorithm] = {'AUC': max_auc, 'fname': file, 'params':params}
        except KeyError:
            best_params[algorithm] = {'AUC': max_auc, 'fname': file, 'params':params}
        
    return best_params


def get_qps_file_for_best_param(filename):
    # take filename and replace 'summary_stats' with 'qps'
    qps_file = filename.replace('summary_stats', 'qps')
    if os.path.exists(qps_file):
        return qps_file
    else:
        return None
def plot_chocolate_side(best_params, qps_files, dataset, output_path):
    plt.figure(figsize=(10, 6))

    for algorithm in best_params:
        qps_file = get_qps_file_for_best_param(best_params[algorithm]['fname'])
        alg_name = best_params[algorithm]['params']


        if not qps_file:
            print(f"No QPS file found for algorithm {algorithm}")
            continue

        qps_df = pd.read_csv(qps_file)
        if 'recall' not in qps_df.columns or 'QPS' not in qps_df.columns:
            print(f"QPS file {qps_file} does not contain 'recall' or 'QPS' columns.")
            continue

        plt.plot(qps_df['recall'], qps_df['QPS'], marker='o', label=f"{alg_name}")

    plt.xlabel('Recall')
    plt.ylabel('QPS')
    plt.title(f'Chocolate Side Plot for Dataset: {dataset}')
    plt.xlim(left=0.8)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Chocolate side plot saved to {output_path}")

def main():
    args = parse_arguments()

    # Flatten the lists (if they are lists of lists)
    summary_stats_files = [item for sublist in args.summary_stats for item in sublist] if any(isinstance(i, list) for i in args.summary_stats) else args.summary_stats
    print("Summary Stats Files:")
    print(summary_stats_files)
    qps_files = [item for sublist in args.qps_files for item in sublist] if any(isinstance(i, list) for i in args.qps_files) else args.qps_files
    print("QPS Files:")
    print(qps_files)

    # Find the best parameters for each algorithm
    best_params = find_best_parameters(summary_stats_files)
    print("Best Parameters:")
    print(best_params)

    # Plot the chocolate side
    plot_chocolate_side(best_params, qps_files, args.dataset, args.output)

if __name__ == "__main__":
    main()
