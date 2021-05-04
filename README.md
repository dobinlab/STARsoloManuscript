# STARsoloManuscript
Code for analyses in the STARsolo manuscript

Directories contain:

data: 10X cell barcode passlists

samples: scRNA-seq data
make -C samples # downloads real data

exe: tool executables
make -C exe     # downloads executables

genomes: genome files and indexes
make -C genomes # downloads genome FASTA/GTF and builds indexes for all tools

sims: simulated data pipeline
make -C sims    # generates simulated data, saves FASTQs in the samples/

count: contains the results of runs
make all        # will run all tools on all datasets. Benchmarking will take a really long time
make real       # will run all tools on the real dataset(s)
make sims       # will run all tolls on the simulated datasets
