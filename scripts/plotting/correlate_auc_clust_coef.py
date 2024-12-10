#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot clustering coefficient vs AUC for each algorithm and parameter combination.")
    parser.add_argument('--dataset', required=True, help='Dataset name.')
    parser.add_argument('--graph_stats', nargs='+', required=True, help='Graph stats CSV files.')
    parser.add_argument('--auc_stats', nargs='+', required=True, help='Summary stats CSV files.')
    parser.add_argument('--output', required=True, help='Output plot file.')
    return parser.parse_args()

def extract_params(filename, dataset):
    basename = os.path.basename(filename)
    # Remove file extension
    remainder = os.path.splitext(basename)[0]

    # Remove known prefixes
    prefixes = ['stats_', 'summary_', f'{dataset}_']
    for prefix in prefixes:
        if remainder.startswith(prefix):
            remainder = remainder[len(prefix):]

    # Remove known suffixes
    suffixes = ['_graph_stats', '_summary', '_graph', '_stats']
    for suffix in suffixes:
        if remainder.endswith(suffix):
            remainder = remainder[:-len(suffix)]

    # Split the remainder into parts
    parts = remainder.split('_')

    algorithms = ['vamana', 'hcnng', 'pyNNDescent']
    algorithm = None
    if parts[0] in algorithms:
        algorithm = parts[0]
        params = parts[1:]
    else:
        # Try to find algorithm in parts
        for i, part in enumerate(parts):
            if part in algorithms:
                algorithm = part
                params = parts[i+1:]
                break
        else:
            print(f"Algorithm not found in filename {filename}")
            return None, None

    param_str = '_'.join(params)
    return algorithm, param_str


def read_and_merge_data(graph_stats_files, summary_stats_files, dataset):
    data = []
    # Create a mapping from (algorithm, param_str) to clustering coefficient
    clust_coef_dict = {}
    for file in graph_stats_files:
        algorithm, param_str = extract_params(file, dataset)
        df = pd.read_csv(file)
        if 'AverageClusteringCoefficient' not in df.columns:
            print(f"Warning: 'AverageClusteringCoefficient' not found in {file}")
            continue
        clust_coef = df['AverageClusteringCoefficient'].iloc[0]
        clust_coef_dict[(algorithm, param_str)] = clust_coef
    # Now, read summary stats and extract AUC values
    for file in summary_stats_files:
        algorithm, param_str = extract_params(file, dataset)
        df = pd.read_csv(file)
        if 'AUC_qps_vs_recall' not in df.columns:
            print(f"Warning: 'AUC_qps_vs_recall' not found in {file}")
            continue
        auc = df['AUC_qps_vs_recall'].iloc[0]
        # Get the clustering coefficient corresponding to this algorithm and param_str
        clust_coef = clust_coef_dict.get((algorithm, param_str), None)
        if clust_coef is not None:
            data.append({
                'Algorithm': algorithm,
                'Parameters': param_str,
                'ClusteringCoefficient': clust_coef,
                'AUC': auc
            })
        else:
            print(f"Warning: No clustering coefficient found for {algorithm} with parameters {param_str}")
    return pd.DataFrame(data)

def plot_clust_coef_vs_auc(df, dataset, output_file):
    plt.figure(figsize=(10, 6))

    algorithms = df['Algorithm'].unique()
    for algorithm in algorithms:
        df_alg = df[df['Algorithm'] == algorithm]
        plt.scatter(df_alg['ClusteringCoefficient'], df_alg['AUC'], label=algorithm)
    
    plt.xlabel('Average Clustering Coefficient')
    plt.ylabel('AUC (QPS vs Recall)')
    plt.title(f'Clustering Coefficient vs AUC for Dataset: {dataset}')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"Plot saved to {output_file}")

def main():
    args = parse_arguments()
    dataset = args.dataset
    graph_stats_files = args.graph_stats
    summary_stats_files = args.auc_stats
    output_file = args.output

    # Read and merge data
    df = read_and_merge_data(graph_stats_files, summary_stats_files, dataset)

    if df.empty:
        print("No data available for plotting.")
        return

    # Generate the scatter plot
    plot_clust_coef_vs_auc(df, dataset, output_file)

if __name__ == "__main__":
    main()
