import argparse
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys


def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate chocolate-side summary for each algorithm.")
    parser.add_argument('--dataset', required=True, help='Dataset name.')
    parser.add_argument('--qps_files', nargs='+', required=True, help='QPS CSV files.')
    parser.add_argument('--output', required=True, help='Output file.')
    parser.add_argument('--k', required=True, type=int, help='Value of k.')
    parser.add_argument('--nthreads', required=True, type=int, help='Number of threads.')
    parser.add_argument('--recalls', nargs='+', type=float, default=[0.8, 0.9, 0.95, 0.99], help='Recall values.')
    return parser.parse_args()


def find_best_parameters(qps_files, k_value, nthreads_value, recalls):
    best_params = {}
    for file in qps_files:
        # Filter files by k and nthreads
        if f"_k{k_value}_" not in file or f"_nthreads{nthreads_value}" not in file:
            continue

        # Extract algorithm name and parameters from filename
        basename = os.path.basename(file)
        dataset = os.path.basename(os.path.dirname(os.path.dirname(file)))
        if basename.startswith(dataset + '_'):
            remainder = basename[len(dataset) + 1:]
        else:
            remainder = basename

        parts = remainder.split('_')
        if len(parts) > 0:
            algorithm = parts[0]
        else:
            print(f"Warning: Unable to extract algorithm name from filename {basename}")
            continue

        # Read the CSV file
        df = pd.read_csv(file)
        if 'recall' not in df.columns or 'QPS' not in df.columns:
            continue

        params = parts[1:-1]  # Exclude algorithm and last part (extension)

        for recall in recalls:
            # Filter data for recall >= threshold
            df_filtered = df[df['recall'] >= recall]
            if df_filtered.empty:
                continue

            # Find the maximum QPS
            # note: this assumes QPS is monotonically decreasing with recall. Fix this. 

            #candQPS should be the first QPS value in the filtered dataframe
            candQPS = df_filtered['QPS'].iloc[0]
            best_rows = df_filtered[df_filtered['QPS'] == candQPS]
            best_row = best_rows.iloc[0]
            qps_at_recall = best_row['QPS']
            recall_at_qps = best_row['recall']

            # Update best_params
            if algorithm not in best_params:
                best_params[algorithm] = {}
            if recall not in best_params[algorithm]:
                best_params[algorithm][recall] = {
                    'QPS': qps_at_recall,
                    'recall': recall_at_qps,
                    'fname': file,
                    'params': params
                }
            else:
                if best_params[algorithm][recall]['QPS'] < qps_at_recall:
                    best_params[algorithm][recall] = {
                        'QPS': qps_at_recall,
                        'recall': recall_at_qps,
                        'fname': file,
                        'params': params
                    }
    return best_params


def write_chocolate_side(best_params, recalls, dataset, output):
    with open(output, 'w') as f:
        f.write(f"# Chocolate-side Summary for Dataset: {dataset}\n\n")
        # Write the table header
        f.write("| Algorithm | " + " | ".join([f"Recall ≥ {recall}" for recall in recalls]) + " |\n")
        f.write("|-----------" + "|------------" * len(recalls) + "|\n")
        for algorithm in best_params:
            row = f"| {algorithm} "
            for recall in recalls:
                if recall in best_params[algorithm]:
                    qps = best_params[algorithm][recall]['QPS']
                    row += f"| {qps:.2f} "
                else:
                    row += "| N/A "
            row += "|\n"
            f.write(row)
        f.write("\n")
        # Optionally, write detailed parameters and files
        f.write("## Detailed Parameters and Files\n\n")
        for algorithm in best_params:
            f.write(f"### {algorithm}\n")
            for recall in recalls:
                if recall in best_params[algorithm]:
                    params = best_params[algorithm][recall]['params']
                    fname = best_params[algorithm][recall]['fname']
                    qps = best_params[algorithm][recall]['QPS']
                    recall_at_qps = best_params[algorithm][recall]['recall']
                    f.write(f"- Recall ≥ {recall}: QPS = {qps:.2f} at recall {recall_at_qps:.2f}\n")
                    f.write(f"  - Params: {' '.join(params)}\n")
                    f.write(f"  - File: {fname}\n")
            f.write("\n")

def main():
    args = parse_arguments()
    k_value = args.k
    nthreads_value = args.nthreads
    recalls = args.recalls

    qps_files = args.qps_files
    if any(isinstance(i, list) for i in qps_files):
        qps_files = [item for sublist in qps_files for item in sublist]

    best_params = find_best_parameters(qps_files, k_value, nthreads_value, recalls)
    write_chocolate_side(best_params, recalls, args.dataset, args.output)

if __name__ == "__main__":
    main()
