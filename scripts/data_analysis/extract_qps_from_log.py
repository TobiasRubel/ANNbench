import csv
import sys
import pandas as pd
from scipy.integrate import simps

def extract_qps_from_log(log_file, csv_file, summary_file):
    with open(log_file, 'r') as infile:
        lines = infile.readlines()

    # Variables to track when the relevant table starts and ends
    start_table = False
    extracted_rows = []
    summary_stats = {}

    # Initialize values for summary statistics
    avg_degree = None
    max_degree = None
    graph_build_time = None

    for line in lines:
        # Capture the average degree and max degree
        if "average degree" in line and "maximum degree" in line:
            parts = line.strip().split(" ")
            avg_degree = float(parts[4])  # Assuming format: "Graph has average degree X and maximum degree Y"
            max_degree = float(parts[-1])

        # Capture graph build time
        if "Graph built in" in line:
            graph_build_time = float(line.strip().split()[3])

        # Start extracting the table after "Graph built in XXX seconds"
        if "Graph built in" in line:
            start_table = True
            continue

        # Stop extracting before "Parlay time: XXX"
        if "Parlay time:" in line:
            break

        # Extract rows only when inside the table
        if start_table and "For" in line:
            parts = line.strip().split(", ")
            row = {}
            for part in parts:
                try:
                    key, value = part.split(" = ")
                except:
                    key, value = part.split(": ")
                if "@" in key:
                    key = key.split(" ")[-1]
                row[key] = float(value)
            extracted_rows.append(row)

    # If no rows were extracted, stop the script
    if not extracted_rows:
        with open(csv_file, 'w', newline='') as csv_out:
            writer = csv.writer(csv_out)
            writer.writerow(["No data extracted from log file"])
        with open(summary_file, 'w', newline='') as summary_out:
            writer = csv.writer(summary_out)
            writer.writerow(["No data extracted from log file"])
        return

    # Writing the extracted table to a CSV file
    with open(csv_file, 'w', newline='') as csv_out:
        fieldnames = extracted_rows[0].keys()
        writer = csv.DictWriter(csv_out, fieldnames=fieldnames)
        
        # Write headers
        writer.writeheader()
        
        # Write data rows
        for row in extracted_rows:
            writer.writerow(row)

    # Calculate AUCs using the trapezoidal rule
    data = pd.DataFrame(extracted_rows)
    auc_qps_vs_recall = simps(data['QPS'], data['recall'])
    auc_visited_vs_recall = simps(data['average visited'], data['recall'])
    auc_cmps_vs_recall = simps(data['average cmps'], data['recall'])

    # Store summary statistics
    summary_stats['average_degree'] = avg_degree
    summary_stats['max_degree'] = max_degree
    summary_stats['graph_build_time'] = graph_build_time
    summary_stats['AUC_qps_vs_recall'] = auc_qps_vs_recall
    summary_stats['AUC_visited_vs_recall'] = auc_visited_vs_recall
    summary_stats['AUC_cmps_vs_recall'] = auc_cmps_vs_recall

    # Write summary statistics to a summary file
    with open(summary_file, 'w', newline='') as summary_out:
        writer = csv.DictWriter(summary_out, fieldnames=summary_stats.keys())
        writer.writeheader()
        writer.writerow(summary_stats)

    print(f"Summary statistics written to {summary_file}")
    print(f"Extracted data written to {csv_file}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_qps_from_log.py <log_file> <csv_file> <summary_file>")
        sys.exit(1)

    log_file = sys.argv[1]
    csv_file = sys.argv[2]
    summary_file = sys.argv[3]

    extract_qps_from_log(log_file, csv_file, summary_file)
