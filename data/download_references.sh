#!/bin/bash

# create conda environment with tools required to build references
baseDir=$(dirname "$0")
cd "$baseDir"
conda env create -f ../task-environment.yml -p conda
source activate ./conda

# download & build references
for REFERENCE in hs37d5viral+GENCODE19 GRCh38viral+GENCODE28 GRCm38viral+GENCODEM25; do
	"$CONDA_PREFIX/var/lib/arriba/download_references.sh" "$REFERENCE"
done

# install database files
mv "$CONDA_PREFIX/var/lib/arriba/"{blacklist,cytobands,known_fusions,protein_domains}* .

# cleanup
rm -rf conda

