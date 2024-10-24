#!/bin/bash

URL = "ftp://ftp.irisa.fr/local/texmex/corpus/sift.tar.gz"

# download data to data directory and decompress

wget $URL -P data/
tar -xvzf data/sift.tar.gz -C data/

# remove tar file
rm data/sift.tar.gz


