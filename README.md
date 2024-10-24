# ANNbench

A collection of scripts for benchmarking Approximate Nearest Neighbors Algorithms on real and synthetic data. 

# Installation & Setup

ParlayANN is a required submodule:

```
git submodule update --init --recursive
```

You'll also need snakemake installed on your system as well as scipy,matplotlib. I'll make the install easier through conda asap. 

Due to the large file sizes, you may prefer to have both the results and data directories be symbolic links to somewhere with a lot of free space. 

Download sift and place it in a data directory following the ParlayANN tutorial for a minimal setup. 

# Running the code

```
snakemake -s main.smk
```

