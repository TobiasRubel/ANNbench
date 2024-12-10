#!/usr/bin/env python3

import argparse
import networkx as nx
import pandas as pd
import os
import struct

def parse_arguments():
    parser = argparse.ArgumentParser(description="Compute clustering coefficient of a graph.")
    parser.add_argument('--input', required=True, help='Input graph file.')
    parser.add_argument('--output', required=True, help='Output stats file.')
    return parser.parse_args()

def read_graph(file_path,directed = False):
    if directed == False:
        G = nx.Graph()
    else:
        G = nx.DiGraph()
    with open(file_path, 'rb') as f:
        # Read the preamble (n and maxDeg)
        n, maxDeg = struct.unpack('II', f.read(8))
        
        # Read sizes array (number of edges for each node)
        sizes = struct.unpack(f'{n}I', f.read(n * 4))
        
        # Read the adjacency list data
        node_index = 0
        for node, size in enumerate(sizes):
            edges = struct.unpack(f'{size}I', f.read(size * 4))
            for edge in edges:
                G.add_edge(node, edge)
                
    return G

def compute_clustering_coefficient(G):
    # Compute the average clustering coefficient
    avg_clustering = nx.approximation.average_clustering(G, trials=1000, seed=0)
    return avg_clustering

def main():
    args = parse_arguments()
    graph_file = args.input
    output_file = args.output

    G = read_graph(graph_file)
    avg_clustering = compute_clustering_coefficient(G)

    # Extract parameters from the filename
    basename = os.path.basename(graph_file)
    dataset, algorithm, *params = basename.replace('.graph', '').split('_')

    # Prepare the data
    data = {
        'Dataset': [dataset],
        'Algorithm': [algorithm],
        'Parameters': ['_'.join(params)],
        'AverageClusteringCoefficient': [avg_clustering]
    }

    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)
    print(f"Clustering coefficient computed and saved to {output_file}")

if __name__ == "__main__":
    main()
