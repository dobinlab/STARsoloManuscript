# requires scanpy installatino
conda activate scanpy
# Preprocess the data for all tools with
python -u real_pbmc_5k.py

# output in write/
# raw -> common -> scale -> neighbors
# neighbors are ready for leiden clustering
