# ANNbench

A collection of scripts for benchmarking Approximate Nearest Neighbors Algorithms on real and synthetic data. 

# Installation & Setup

ParlayANN is a required submodule:

```
git submodule init
git submodule update
```

You'll also need snakemake installed on your system. 

Due to the large file sizes, you may prefer to have both the results and data directories be symbolic links to somewhere with a lot of free space. 
