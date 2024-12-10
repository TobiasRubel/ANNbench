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
    parser.add_argument('--k', required=True, type=int, help='Value of k.')
    parser.add_argument('--nthreads', required=True, type=int, help='Number of threads.')
    return parser.parse_args()


def find_best_parameters(summary_stats_files, k_value, nthreads_value):
    AUC = 'AUC_qps_vs_recall'
    best_params = {}
    for file in summary_stats_files:
        # Filter files by k and nthreads
        if f"_k{k_value}_" not in file or f"_nthreads{nthreads_value}" not in file:
            continue
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
        params = parts[0:-2]

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
        
# def plot_chocolate_side(best_params, qps_files, dataset, output_path):
#     plt.figure(figsize=(10, 6))

#     for algorithm in best_params:
#         qps_file = get_qps_file_for_best_param(best_params[algorithm]['fname'])
#         alg_name = best_params[algorithm]['params']


#         if not qps_file:
#             print(f"No QPS file found for algorithm {algorithm}")
#             continue

#         qps_df = pd.read_csv(qps_file)
#         if 'recall' not in qps_df.columns or 'QPS' not in qps_df.columns:
#             print(f"QPS file {qps_file} does not contain 'recall' or 'QPS' columns.")
#             continue

#         plt.plot(qps_df['recall'], qps_df['QPS'], marker='o', label=f"{alg_name}")

#     plt.xlabel('Recall')
#     plt.ylabel('QPS')
#     plt.title(f'Chocolate Side Plot for Dataset: {dataset}')
#     minx = 0.8
#     plt.xlim(left=minx)
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig(output_path)
#     plt.close()
#     print(f"Chocolate side plot saved to {output_path}")

def plot_chocolate_side(best_params, qps_files, dataset, output_path):
    plt.figure(figsize=(10, 6))
    minx = 0.8
    max_qps = 0  # To keep track of the maximum QPS across all algorithms

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

        # Filter data to recall ≥ minx
        qps_df_filtered = qps_df[qps_df['recall'] >= minx]
        if qps_df_filtered.empty:
            print(f"No data for algorithm {algorithm} with recall ≥ {minx}")
            continue

        # Update max_qps with the maximum QPS from the filtered data
        current_max_qps = qps_df_filtered['QPS'].max()
        if current_max_qps > max_qps:
            max_qps = current_max_qps

        # Plot the filtered data
        plt.plot(qps_df_filtered['recall'], qps_df_filtered['QPS'], marker='o', label=f"{alg_name}")

    # Check if there is data to plot
    if max_qps == 0:
        print(f"No data to plot for recall ≥ {minx}")
        return

    plt.xlabel('Recall')
    plt.ylabel('QPS')
    plt.title(f'Chocolate Side Plot for Dataset: {dataset}')
    plt.xlim(left=minx)
    # Set y-axis limits with a small margin above max_qps
    plt.ylim(bottom=0, top=max_qps * 1.05)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Chocolate side plot saved to {output_path}")

def main():
    args = parse_arguments()
    k_value = args.k
    nthreads_value = args.nthreads

    # Flatten the lists (if they are lists of lists)
    summary_stats_files = [item for sublist in args.summary_stats for item in sublist] if any(isinstance(i, list) for i in args.summary_stats) else args.summary_stats
    qps_files = [item for sublist in args.qps_files for item in sublist] if any(isinstance(i, list) for i in args.qps_files) else args.qps_files

    # Filter qps_files by k and nthreads
    qps_files = [f for f in qps_files if f"_k{k_value}_" in f and f"_nthreads{nthreads_value}" in f]

    # Find the best parameters for each algorithm
    best_params = find_best_parameters(summary_stats_files, k_value, nthreads_value)

    # Plot the chocolate side
    plot_chocolate_side(best_params, qps_files, args.dataset, args.output)

if __name__ == "__main__":
    main()
