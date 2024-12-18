import os
import yaml
import glob

# Load configuration
configfile: "config.yaml"

# Directories
PARLAY_DIR = "ParlayANN/algorithms"
RESULTS_DIR = "results"
DATA_DIR = "data"

# Dataset names and global parameters
datasets = config["datasets"].keys()
k_values = [10]  # Change this for param sweeps
nthreads_values = [32]  # Change this for param sweeps


algorithms = ["hcnng", "vamana", "pyNNDescent"]

# Parameters for HCNNG
HCNNG_cluster_sizes = [500, 1000, 2000]
#HCNNG_cluster_sizes = [1000]
HCNNG_num_clusters_values = [20,30,40]
#HCNNG_num_clusters_values = [30]
# Parameters for Vamana
vamana_R_values = [16,32,48]
#vamana_R_values = [32]
vamana_L_values = [32,64,128]
#vamana_L_values = [64]
vamana_alpha_values = [1.0, 1.2, 1.4]
#vamana_alpha_values = [1.2]
# Parameters for HNSW

HNSW_M_values = [20]
HNSW_efc_values = [50]
HNSW_alpha_values = [0.9]
HNSW_ml_values = [34]

# Parameters for pyNNDescent

#pyNNDescent_R_values = [20, 30]
pyNNDescent_R_values = [30]
pyNNDescent_cluster_sizes = [100]
#pyNNDescent_num_clusters_values = [5,10,20]
pyNNDescent_num_clusters_values = [10]
#pyNNDescent_alpha_values = [1.2,1.4]
pyNNDescent_alpha_values = [1.2]
#pyNNDescent_delta_values = [0.02, 0.05]
pyNNDescent_delta_values = [0.02]


# Rule for the final outputs
rule all:
    input:
        # create targets for .log files from running algorithms
        expand(os.path.join(RESULTS_DIR, "{dataset}","logs", "{dataset}_hcnng_cluster{cluster}_num_clusters{num_clusters}_k{k}_nthreads{nthreads}.log"),
               cluster=HCNNG_cluster_sizes, num_clusters=HCNNG_num_clusters_values, k=k_values, nthreads=nthreads_values, dataset=datasets),
        expand(os.path.join(RESULTS_DIR, "{dataset}","logs", "{dataset}_vamana_R{R}_L{L}_alpha{alpha}_k{k}_nthreads{nthreads}.log"),
               R=vamana_R_values, L=vamana_L_values, alpha=vamana_alpha_values, k=k_values, nthreads=nthreads_values, dataset=datasets),
        expand(os.path.join(RESULTS_DIR, "{dataset}","logs", "{dataset}_pyNNDescent_R{R}_cluster{cluster}_num_clusters{num_clusters}_alpha{alpha}_delta{delta}_k{k}_nthreads{nthreads}.log"),
               R=pyNNDescent_R_values, cluster=pyNNDescent_cluster_sizes, num_clusters=pyNNDescent_num_clusters_values, alpha=pyNNDescent_alpha_values, delta=pyNNDescent_delta_values, k=k_values, nthreads=nthreads_values, dataset=datasets),
        
        #create targets for graph files
        expand(os.path.join(DATA_DIR, "{dataset}","graphs", "{dataset}_hcnng_cluster{cluster}_num_clusters{num_clusters}_k{k}_nthreads{nthreads}.graph"),
               cluster=HCNNG_cluster_sizes, num_clusters=HCNNG_num_clusters_values, k=k_values, nthreads=nthreads_values, dataset=datasets),
        expand(os.path.join(DATA_DIR, "{dataset}","graphs", "{dataset}_vamana_R{R}_L{L}_alpha{alpha}_k{k}_nthreads{nthreads}.graph"),
               R=vamana_R_values, L=vamana_L_values, alpha=vamana_alpha_values, k=k_values, nthreads=nthreads_values, dataset=datasets),
        expand(os.path.join(DATA_DIR, "{dataset}","graphs", "{dataset}_pyNNDescent_R{R}_cluster{cluster}_num_clusters{num_clusters}_alpha{alpha}_delta{delta}_k{k}_nthreads{nthreads}.graph"),
               R=pyNNDescent_R_values, cluster=pyNNDescent_cluster_sizes, num_clusters=pyNNDescent_num_clusters_values, alpha=pyNNDescent_alpha_values, delta=pyNNDescent_delta_values, k=k_values, nthreads=nthreads_values, dataset=datasets),



        # create targets for .csv files from converting .log files into more parsable format
        expand(os.path.join(RESULTS_DIR, "{dataset}","qps", "{dataset}_hcnng_cluster{cluster}_num_clusters{num_clusters}_k{k}_nthreads{nthreads}.csv"),
               cluster=HCNNG_cluster_sizes, num_clusters=HCNNG_num_clusters_values, k=k_values, nthreads=nthreads_values, dataset=datasets),
        expand(os.path.join(RESULTS_DIR, "{dataset}","qps", "{dataset}_vamana_R{R}_L{L}_alpha{alpha}_k{k}_nthreads{nthreads}.csv"),
               R=vamana_R_values, L=vamana_L_values, alpha=vamana_alpha_values, k=k_values, nthreads=nthreads_values, dataset=datasets),
        expand(os.path.join(RESULTS_DIR, "{dataset}","qps", "{dataset}_pyNNDescent_R{R}_cluster{cluster}_num_clusters{num_clusters}_alpha{alpha}_delta{delta}_k{k}_nthreads{nthreads}.csv"),
               R=pyNNDescent_R_values, cluster=pyNNDescent_cluster_sizes, num_clusters=pyNNDescent_num_clusters_values, alpha=pyNNDescent_alpha_values, delta=pyNNDescent_delta_values, k=k_values, nthreads=nthreads_values, dataset=datasets),

        # create targets for generating summary statistics
        expand(os.path.join(RESULTS_DIR, "{dataset}","summary_stats", "{dataset}_hcnng_cluster{cluster}_num_clusters{num_clusters}_k{k}_nthreads{nthreads}.csv"),
               cluster=HCNNG_cluster_sizes, num_clusters=HCNNG_num_clusters_values, k=k_values, nthreads=nthreads_values, dataset=datasets),
        expand(os.path.join(RESULTS_DIR, "{dataset}","summary_stats", "{dataset}_vamana_R{R}_L{L}_alpha{alpha}_k{k}_nthreads{nthreads}.csv"),
               R=vamana_R_values, L=vamana_L_values, alpha=vamana_alpha_values, k=k_values, nthreads=nthreads_values, dataset=datasets),
        expand(os.path.join(RESULTS_DIR, "{dataset}","summary_stats", "{dataset}_pyNNDescent_R{R}_cluster{cluster}_num_clusters{num_clusters}_alpha{alpha}_delta{delta}_k{k}_nthreads{nthreads}.csv"),
               R=pyNNDescent_R_values, cluster=pyNNDescent_cluster_sizes, num_clusters=pyNNDescent_num_clusters_values, alpha=pyNNDescent_alpha_values, delta=pyNNDescent_delta_values, k=k_values, nthreads=nthreads_values, dataset=datasets),

        # Chocolate-side plots
        expand(os.path.join(RESULTS_DIR, "{dataset}", "plots", "{dataset}_k{k}_nthreads{nthreads}_chocolate_side.png"),
               dataset=datasets, k=k_values, nthreads=nthreads_values),
        # Chocolate summary statistics
        expand(os.path.join(RESULTS_DIR, "{dataset}", "summary_stats", "{dataset}_k{k}_nthreads{nthreads}_chocolate_summary.md"),
               dataset=datasets, k=k_values, nthreads=nthreads_values),
        
        # clustering coefficient extimations
        expand(os.path.join(RESULTS_DIR, "{dataset}", "summary_stats", "clustering_coefficient.txt"),
               dataset=datasets),
        expand(os.path.join(RESULTS_DIR, "{dataset}", "plots", "clustering_coefficient.png"),
               dataset=datasets)

        # # HCNNG graph stats
        #     expand(
        #     os.path.join(RESULTS_DIR, "{dataset}", "summary_stats",
        #                  "graph_stats_{dataset}_hcnng_cluster{cluster}_num_clusters{num_clusters}_k{k}_nthreads{nthreads}.csv"),
        #     dataset=datasets,
        #     cluster=HCNNG_cluster_sizes,
        #     num_clusters=HCNNG_num_clusters_values,
        #     k=k_values,
        #     nthreads=nthreads_values
        # ),
        # # Vamana graph stats
        # expand(
        #     os.path.join(RESULTS_DIR, "{dataset}", "summary_stats",
        #                  "graph_stats_{dataset}_vamana_R{R}_L{L}_alpha{alpha}_k{k}_nthreads{nthreads}.csv"),
        #     dataset=datasets,
        #     R=vamana_R_values,
        #     L=vamana_L_values,
        #     alpha=vamana_alpha_values,
        #     k=k_values,
        #     nthreads=nthreads_values
        # ),
        # # pyNNDescent graph stats
        # expand(
        #     os.path.join(RESULTS_DIR, "{dataset}", "summary_stats",
        #                  "graph_stats_{dataset}_pyNNDescent_R{R}_cluster{cluster}_num_clusters{num_clusters}_alpha{alpha}_delta{delta}_k{k}_nthreads{nthreads}.csv"),
        #     dataset=datasets,
        #     R=pyNNDescent_R_values,
        #     cluster=pyNNDescent_cluster_sizes,
        #     num_clusters=pyNNDescent_num_clusters_values,
        #     alpha=pyNNDescent_alpha_values,
        #     delta=pyNNDescent_delta_values,
        #     k=k_values,
        #     nthreads=nthreads_values
        # ),

        # # Correlation between AUC and clustering coefficient
        # expand(os.path.join(RESULTS_DIR, "{dataset}", "plots", "{dataset}_k{k}_nthreads{nthreads}_auc_vs_clust.png"),
        #        dataset=datasets, k=k_values, nthreads=nthreads_values),

# rules for building code
        
rule build_hcnng:
    output:
        os.path.join(PARLAY_DIR, "HCNNG", "neighbors")
    shell:
        "cd {PARLAY_DIR}/HCNNG && make"

rule build_vamana:
    output:
        os.path.join(PARLAY_DIR, "vamana", "neighbors")
    shell:
        "cd {PARLAY_DIR}/vamana && make"

rule build_pynndescent:
    output:
        os.path.join(PARLAY_DIR, "pyNNDescent", "neighbors")
    shell:
        "cd {PARLAY_DIR}/pyNNDescent && make"

rule build_clustering_coefficient:
    output:
        os.path.join("scripts", "data_analysis", "compute_clustering_coeffs")
    shell:
        "cd scripts/data_analysis && make"

# rules for running programs

rule run_hcnng:
    input:
        neighbors=os.path.join(PARLAY_DIR, "HCNNG", "neighbors"),
        query_file=lambda wildcards: os.path.join(DATA_DIR, wildcards.dataset, config["datasets"][wildcards.dataset]["query_file"]),
        base_file=lambda wildcards: os.path.join(DATA_DIR, wildcards.dataset, config["datasets"][wildcards.dataset]["base_file"]),
        gt_path=lambda wildcards: os.path.join(DATA_DIR, wildcards.dataset, config["datasets"][wildcards.dataset]["gt_file"])
    output:
        log_path=os.path.join(RESULTS_DIR, "{dataset}", "logs", "{dataset}_hcnng_cluster{cluster}_num_clusters{num_clusters}_k{k}_nthreads{nthreads}.log"),
        graph_path=os.path.join(DATA_DIR, "{dataset}", "graphs", "{dataset}_hcnng_cluster{cluster}_num_clusters{num_clusters}_k{k}_nthreads{nthreads}.graph")
    params:
        data_type=lambda wildcards: config["datasets"][wildcards.dataset]["data_type"],
        dist_func=lambda wildcards: config["datasets"][wildcards.dataset]["dist_func"],
    shell:
        """
        PARLAY_NUM_THREADS={wildcards.nthreads} {input.neighbors} \
        -cluster_size {wildcards.cluster} -mst_deg 3 -num_clusters {wildcards.num_clusters} \
        -graph_outfile {output.graph_path} -query_path {input.query_file} -gt_path {input.gt_path} \
        -data_type {params.data_type} -dist_func {params.dist_func} -base_path {input.base_file} -k {wildcards.k} \
        > {output.log_path} 2>&1
        """

rule run_vamana:
    input:
        neighbors=os.path.join(PARLAY_DIR, "vamana", "neighbors"),
        query_file=lambda wildcards: os.path.join(DATA_DIR, wildcards.dataset, config["datasets"][wildcards.dataset]["query_file"]),
        base_file=lambda wildcards: os.path.join(DATA_DIR, wildcards.dataset, config["datasets"][wildcards.dataset]["base_file"]),
        gt_path=lambda wildcards: os.path.join(DATA_DIR, wildcards.dataset, config["datasets"][wildcards.dataset]["gt_file"])
    output:
        log_path=os.path.join(RESULTS_DIR, "{dataset}", "logs", "{dataset}_vamana_R{R}_L{L}_alpha{alpha}_k{k}_nthreads{nthreads}.log"),
        graph_path=os.path.join(DATA_DIR, "{dataset}", "graphs", "{dataset}_vamana_R{R}_L{L}_alpha{alpha}_k{k}_nthreads{nthreads}.graph")
    params:
        data_type=lambda wildcards: config["datasets"][wildcards.dataset]["data_type"],
        dist_func=lambda wildcards: config["datasets"][wildcards.dataset]["dist_func"],
    shell:
        """
        PARLAY_NUM_THREADS={wildcards.nthreads} {input.neighbors} \
        -R {wildcards.R} -L {wildcards.L} -alpha {wildcards.alpha} \
        -graph_outfile {output.graph_path} -query_path {input.query_file} -gt_path {input.gt_path} \
        -data_type {params.data_type} -dist_func {params.dist_func} -base_path {input.base_file} -k {wildcards.k} \
        > {output.log_path} 2>&1
        """

rule run_pynndescent:
    input:
        neighbors=os.path.join(PARLAY_DIR, "pyNNDescent", "neighbors"),
        query_file=lambda wildcards: os.path.join(DATA_DIR, wildcards.dataset, config["datasets"][wildcards.dataset]["query_file"]),
        base_file=lambda wildcards: os.path.join(DATA_DIR, wildcards.dataset, config["datasets"][wildcards.dataset]["base_file"]),
        gt_path=lambda wildcards: os.path.join(DATA_DIR, wildcards.dataset, config["datasets"][wildcards.dataset]["gt_file"])
    output:
        log_path=os.path.join(RESULTS_DIR, "{dataset}", "logs", "{dataset}_pyNNDescent_R{R}_cluster{cluster}_num_clusters{num_clusters}_alpha{alpha}_delta{delta}_k{k}_nthreads{nthreads}.log"),
        graph_path=os.path.join(DATA_DIR, "{dataset}", "graphs", "{dataset}_pyNNDescent_R{R}_cluster{cluster}_num_clusters{num_clusters}_alpha{alpha}_delta{delta}_k{k}_nthreads{nthreads}.graph")
    params:
        data_type=lambda wildcards: config["datasets"][wildcards.dataset]["data_type"],
        dist_func=lambda wildcards: config["datasets"][wildcards.dataset]["dist_func"],
    shell:
        """
        PARLAY_NUM_THREADS={wildcards.nthreads} {input.neighbors} \
        -R {wildcards.R} -cluster_size {wildcards.cluster} -num_clusters {wildcards.num_clusters} -alpha {wildcards.alpha} -delta {wildcards.delta} \
        -graph_outfile {output.graph_path} -query_path {input.query_file} -gt_path {input.gt_path} \
        -data_type {params.data_type} -dist_func {params.dist_func} -base_path {input.base_file} -k {wildcards.k} \
        > {output.log_path} 2>&1
        """

# rules for converting logs to csv

rule hcnng_log_to_csv:
    input:
        log_file=os.path.join(RESULTS_DIR, "{dataset}/logs/", "{dataset}_hcnng_cluster{cluster}_num_clusters{num_clusters}_k{k}_nthreads{nthreads}.log")
    output:
        csv_file=os.path.join(RESULTS_DIR, "{dataset}/qps/", "{dataset}_hcnng_cluster{cluster}_num_clusters{num_clusters}_k{k}_nthreads{nthreads}.csv"),
        summary_file=os.path.join(RESULTS_DIR, "{dataset}/summary_stats/", "{dataset}_hcnng_cluster{cluster}_num_clusters{num_clusters}_k{k}_nthreads{nthreads}.csv")
    shell:
        """
        python3 scripts/data_analysis/extract_qps_from_log.py {input.log_file} {output.csv_file} {output.summary_file}
        """

rule vamana_log_to_csv:
    input:
        log_file=os.path.join(RESULTS_DIR, "{dataset}/logs/", "{dataset}_vamana_R{R}_L{L}_alpha{alpha}_k{k}_nthreads{nthreads}.log")
    output:
        csv_file=os.path.join(RESULTS_DIR, "{dataset}/qps/", "{dataset}_vamana_R{R}_L{L}_alpha{alpha}_k{k}_nthreads{nthreads}.csv"),
        summary_file=os.path.join(RESULTS_DIR, "{dataset}/summary_stats/", "{dataset}_vamana_R{R}_L{L}_alpha{alpha}_k{k}_nthreads{nthreads}.csv")
    shell:
        """
        python3 scripts/data_analysis/extract_qps_from_log.py {input.log_file} {output.csv_file} {output.summary_file}
        """

rule pynndescent_log_to_csv:
    input:
        log_file=os.path.join(RESULTS_DIR, "{dataset}/logs/", "{dataset}_pyNNDescent_R{R}_cluster{cluster}_num_clusters{num_clusters}_alpha{alpha}_delta{delta}_k{k}_nthreads{nthreads}.log")
    output:
        csv_file=os.path.join(RESULTS_DIR, "{dataset}/qps/", "{dataset}_pyNNDescent_R{R}_cluster{cluster}_num_clusters{num_clusters}_alpha{alpha}_delta{delta}_k{k}_nthreads{nthreads}.csv"),
        summary_file=os.path.join(RESULTS_DIR, "{dataset}/summary_stats/", "{dataset}_pyNNDescent_R{R}_cluster{cluster}_num_clusters{num_clusters}_alpha{alpha}_delta{delta}_k{k}_nthreads{nthreads}.csv"),
    shell:
        """
        python3 scripts/data_analysis/extract_qps_from_log.py {input.log_file} {output.csv_file} {output.summary_file}
        """

# generate graph summary statistics

# rule compute_graph_stats_hcnng:
#     input:
#         graph_file=os.path.join(DATA_DIR, "{dataset}/graphs/", "{dataset}_hcnng_cluster{cluster}_num_clusters{num_clusters}_k{k}_nthreads{nthreads}.graph")
#     output:
#         os.path.join(
#             RESULTS_DIR, "{dataset}", "summary_stats",
#             "graph_stats_{dataset}_hcnng_cluster{cluster}_num_clusters{num_clusters}_k{k}_nthreads{nthreads}.csv"
#         )
#     params:
#         script="scripts/data_analysis/compute_clustering_coefficient.py"
#     shell:
#         """
#         python {params.script} --input {input.graph_file} --output {output}
#         """


# rule compute_graph_stats_vamana:
#     input:
#         graph_file=os.path.join(DATA_DIR, "{dataset}/graphs/", "{dataset}_vamana_R{R}_L{L}_alpha{alpha}_k{k}_nthreads{nthreads}.graph")

#     output:
#         os.path.join(
#             RESULTS_DIR, "{dataset}", "summary_stats",
#             "graph_stats_{dataset}_vamana_R{R}_L{L}_alpha{alpha}_k{k}_nthreads{nthreads}.csv"
#         )
#     params:
#         script="scripts/data_analysis/compute_clustering_coefficient.py"
#     shell:
#         """
#         python {params.script} --input {input.graph_file} --output {output}
#         """


# rule compute_graph_stats_pynndescent:
#     input:
#         graph_file=os.path.join(DATA_DIR, "{dataset}/graphs/", "{dataset}_pyNNDescent_R{R}_cluster{cluster}_num_clusters{num_clusters}_alpha{alpha}_delta{delta}_k{k}_nthreads{nthreads}.graph")
#     output:
#         os.path.join(
#             RESULTS_DIR, "{dataset}", "summary_stats",
#             "graph_stats_{dataset}_pyNNDescent_R{R}_cluster{cluster}_num_clusters{num_clusters}_alpha{alpha}_delta{delta}_k{k}_nthreads{nthreads}.csv"
#         )
#     params:
#         script="scripts/data_analysis/compute_clustering_coefficient.py"
#     shell:
#         """
#         python {params.script} --input {input.graph_file} --output {output}
#         """




# rules for plotting


# Rule to plot dataset QPS from CSV files
# rule plot_dataset_qps:
#     input:
#         csv_files=lambda wildcards: glob.glob(os.path.join(RESULTS_DIR, wildcards.dataset, "qps", f"*k{wildcards.k}_*.csv"))
#     output:
#         plot_file=os.path.join(RESULTS_DIR, "{dataset}/plots/", "{dataset}_k{k}_nthreads{nthreads}_qps.pdf")
#     shell:
#         """
#         python3 -W ignore scripts/plotting/plot_qps.py {input.csv_files} {output.plot_file}
#         """


# Rule to plot the "chocolate side"
# which here is defined as the QPS vs recall for the highest AUC parameter combination for an algorithm.
def flatten(list_of_lists):
    return [item for sublist in list_of_lists for item in sublist]

rule plot_chocolate_side:
    input:
        summary_stats=lambda wildcards: flatten([
            # HCNNG
            expand(
                os.path.join(RESULTS_DIR, wildcards.dataset, "summary_stats", "{dataset}_hcnng_cluster{cluster}_num_clusters{num_clusters}_k{k}_nthreads{nthreads}.csv"),
                cluster=HCNNG_cluster_sizes, num_clusters=HCNNG_num_clusters_values,
                k=[wildcards.k], nthreads=[wildcards.nthreads], dataset=wildcards.dataset
            ),
            # Vamana
            expand(
                os.path.join(RESULTS_DIR, wildcards.dataset, "summary_stats", "{dataset}_vamana_R{R}_L{L}_alpha{alpha}_k{k}_nthreads{nthreads}.csv"),
                R=vamana_R_values, L=vamana_L_values, alpha=vamana_alpha_values,
                k=[wildcards.k], nthreads=[wildcards.nthreads], dataset=wildcards.dataset
            ),
            # pyNNDescent
            expand(
                os.path.join(RESULTS_DIR, wildcards.dataset, "summary_stats", "{dataset}_pyNNDescent_R{R}_cluster{cluster}_num_clusters{num_clusters}_alpha{alpha}_delta{delta}_k{k}_nthreads{nthreads}.csv"),
                R=pyNNDescent_R_values, cluster=pyNNDescent_cluster_sizes,
                num_clusters=pyNNDescent_num_clusters_values, alpha=pyNNDescent_alpha_values,
                delta=pyNNDescent_delta_values, k=[wildcards.k], nthreads=[wildcards.nthreads], dataset=wildcards.dataset
            )
        ]),
        qps_files=lambda wildcards: flatten([
            # HCNNG
            expand(
                os.path.join(RESULTS_DIR, wildcards.dataset, "qps", "{dataset}_hcnng_cluster{cluster}_num_clusters{num_clusters}_k{k}_nthreads{nthreads}.csv"),
                cluster=HCNNG_cluster_sizes, num_clusters=HCNNG_num_clusters_values,
                k=[wildcards.k], nthreads=[wildcards.nthreads], dataset=wildcards.dataset
            ),
            # Vamana
            expand(
                os.path.join(RESULTS_DIR, wildcards.dataset, "qps", "{dataset}_vamana_R{R}_L{L}_alpha{alpha}_k{k}_nthreads{nthreads}.csv"),
                R=vamana_R_values, L=vamana_L_values, alpha=vamana_alpha_values,
                k=[wildcards.k], nthreads=[wildcards.nthreads], dataset=wildcards.dataset
            ),
            # pyNNDescent
            expand(
                os.path.join(RESULTS_DIR, wildcards.dataset, "qps", "{dataset}_pyNNDescent_R{R}_cluster{cluster}_num_clusters{num_clusters}_alpha{alpha}_delta{delta}_k{k}_nthreads{nthreads}.csv"),
                R=pyNNDescent_R_values, cluster=pyNNDescent_cluster_sizes,
                num_clusters=pyNNDescent_num_clusters_values, alpha=pyNNDescent_alpha_values,
                delta=pyNNDescent_delta_values, k=[wildcards.k], nthreads=[wildcards.nthreads], dataset=wildcards.dataset
            )
        ]),
    output:
        os.path.join(RESULTS_DIR, "{dataset}", "plots", "{dataset}_k{k}_nthreads{nthreads}_chocolate_side.png")
    params:
        plot_script="scripts/plotting/plot_chocolate_side.py"
    shell:
        """
        python {params.plot_script} --dataset {wildcards.dataset} \
            --k {wildcards.k} --nthreads {wildcards.nthreads} \
            --summary_stats {input.summary_stats} \
            --qps_files {input.qps_files} \
            --output {output}
        """


# Chocolate Summary Statistics
# Where we just take the highest QPS at a given recall for each algorithm. Here we'll just use 0.9 recall for simplicity at the moment.
summary_recall_values = [0.8, 0.9, 0.95, 0.99]

rule chocolate_summary_stats:
    input:
        qps_files=lambda wildcards: flatten([
            # HCNNG
            expand(
                os.path.join(RESULTS_DIR, wildcards.dataset, "qps", "{dataset}_hcnng_cluster{cluster}_num_clusters{num_clusters}_k{k}_nthreads{nthreads}.csv"),
                cluster=HCNNG_cluster_sizes, num_clusters=HCNNG_num_clusters_values,
                k=[wildcards.k], nthreads=[wildcards.nthreads], dataset=wildcards.dataset
            ),
            # Vamana
            expand(
                os.path.join(RESULTS_DIR, wildcards.dataset, "qps", "{dataset}_vamana_R{R}_L{L}_alpha{alpha}_k{k}_nthreads{nthreads}.csv"),
                R=vamana_R_values, L=vamana_L_values, alpha=vamana_alpha_values,
                k=[wildcards.k], nthreads=[wildcards.nthreads], dataset=wildcards.dataset
            ),
            # pyNNDescent
            expand(
                os.path.join(RESULTS_DIR, wildcards.dataset, "qps", "{dataset}_pyNNDescent_R{R}_cluster{cluster}_num_clusters{num_clusters}_alpha{alpha}_delta{delta}_k{k}_nthreads{nthreads}.csv"),
                R=pyNNDescent_R_values, cluster=pyNNDescent_cluster_sizes,
                num_clusters=pyNNDescent_num_clusters_values, alpha=pyNNDescent_alpha_values,
                delta=pyNNDescent_delta_values, k=[wildcards.k], nthreads=[wildcards.nthreads], dataset=wildcards.dataset
            )
        ]),


    output:
        os.path.join(RESULTS_DIR, "{dataset}", "summary_stats", "{dataset}_k{k}_nthreads{nthreads}_chocolate_summary.md")
    params:
        summary_script="scripts/data_analysis/chocolate_summary.py",
        recalls=' '.join(map(str, summary_recall_values))
    shell:
        """
        python {params.summary_script} --dataset {wildcards.dataset} \
            --k {wildcards.k} --nthreads {wildcards.nthreads} \
            --qps_files {input.qps_files} \
            --recalls {params.recalls} \
            --output {output}
        """


rule correlate_auc_clust_coef:
    input:
        auc_stats=lambda wildcards: flatten([
            # HCNNG
            expand(
                os.path.join(RESULTS_DIR, wildcards.dataset, "summary_stats", "{dataset}_hcnng_cluster{cluster}_num_clusters{num_clusters}_k{k}_nthreads{nthreads}.csv"),
                cluster=HCNNG_cluster_sizes, num_clusters=HCNNG_num_clusters_values,
                k=[wildcards.k], nthreads=[wildcards.nthreads], dataset=wildcards.dataset
            ),
            # Vamana
            expand(
                os.path.join(RESULTS_DIR, wildcards.dataset, "summary_stats", "{dataset}_vamana_R{R}_L{L}_alpha{alpha}_k{k}_nthreads{nthreads}.csv"),
                R=vamana_R_values, L=vamana_L_values, alpha=vamana_alpha_values,
                k=[wildcards.k], nthreads=[wildcards.nthreads], dataset=wildcards.dataset
            ),
            # pyNNDescent
            expand(
                os.path.join(RESULTS_DIR, wildcards.dataset, "summary_stats", "{dataset}_pyNNDescent_R{R}_cluster{cluster}_num_clusters{num_clusters}_alpha{alpha}_delta{delta}_k{k}_nthreads{nthreads}.csv"),
                R=pyNNDescent_R_values, cluster=pyNNDescent_cluster_sizes,
                num_clusters=pyNNDescent_num_clusters_values, alpha=pyNNDescent_alpha_values,
                delta=pyNNDescent_delta_values, k=[wildcards.k], nthreads=[wildcards.nthreads], dataset=wildcards.dataset
            )
        ]),
        graph_stats=lambda wildcards: flatten([
            # HCNNG
            expand(
                os.path.join(RESULTS_DIR, wildcards.dataset, "summary_stats", "graph_stats_{dataset}_hcnng_cluster{cluster}_num_clusters{num_clusters}_k{k}_nthreads{nthreads}.csv"),
                cluster=HCNNG_cluster_sizes, num_clusters=HCNNG_num_clusters_values,
                k=[wildcards.k], nthreads=[wildcards.nthreads], dataset=wildcards.dataset
            ),
            # Vamana
            expand(
                os.path.join(RESULTS_DIR, wildcards.dataset, "summary_stats", "graph_stats_{dataset}_vamana_R{R}_L{L}_alpha{alpha}_k{k}_nthreads{nthreads}.csv"),
                R=vamana_R_values, L=vamana_L_values, alpha=vamana_alpha_values,
                k=[wildcards.k], nthreads=[wildcards.nthreads], dataset=wildcards.dataset
            ),
            # pyNNDescent
            expand(
                os.path.join(RESULTS_DIR, wildcards.dataset, "summary_stats", "graph_stats_{dataset}_pyNNDescent_R{R}_cluster{cluster}_num_clusters{num_clusters}_alpha{alpha}_delta{delta}_k{k}_nthreads{nthreads}.csv"),
                R=pyNNDescent_R_values, cluster=pyNNDescent_cluster_sizes,
                num_clusters=pyNNDescent_num_clusters_values, alpha=pyNNDescent_alpha_values,
                delta=pyNNDescent_delta_values, k=[wildcards.k], nthreads=[wildcards.nthreads], dataset=wildcards.dataset
            )
        ]),
    output:
        os.path.join(RESULTS_DIR, "{dataset}", "plots", "{dataset}_k{k}_nthreads{nthreads}_auc_vs_clust.png")
    params:
        plot_script="scripts/plotting/correlate_auc_clust_coef.py"
    shell:
        """
        python {params.plot_script} --dataset {wildcards.dataset} \
            --auc_stats {input.auc_stats} \
            --graph_stats {input.graph_stats} \
            --output {output}
        """


# rules for generating clustering coefficient stuffs

rule compute_clustering_coefficient:
    input:
        base_dir=lambda wildcards: os.path.join(DATA_DIR, wildcards.dataset, config["datasets"][wildcards.dataset]["base_file"])
    output:
        os.path.join(RESULTS_DIR, "{dataset}", "summary_stats", "clustering_coefficient.txt")
    params:
        executable=os.path.join("scripts", "data_analysis", "compute_clustering_coeffs")
    shell:
        """
        PARLAY_NUM_THREADS=32 {params.executable} {input.base_dir} {output}
        """


rule plot_clustering_coefficient:
    input:
        os.path.join(RESULTS_DIR, "{dataset}", "summary_stats", "clustering_coefficient.txt")
    output:
        os.path.join(RESULTS_DIR, "{dataset}", "plots", "clustering_coefficient.png")
    params:
        executable=os.path.join("scripts", "plotting", "plot_clustering_coeffs.py")
    shell:
        """
        python {params.executable} {input} {output}
        """
