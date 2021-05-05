#%% adapted from Dinar Yunusov's code: 03/22/2021

#%% imports
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import os # to run bash commands from python
import math # to calculate log2

from sklearn.cluster import KMeans
import sklearn.metrics as sm

# new:
# remove regressions
# hvg: flavor='seurat_v3', 3000 top genes
# use neighbors_dense for reproducible NN calculation: not affected by random_state

#%% variables
mtxPrefix = '/scratch/dobin/STAR/STARsoloPreprint/count/'
genome1 = 'mouse_CR_3.0.0'
sample1 = '10X/3/pbmc_5k/20/run1'
write_dir = './write_' + 'pbmc5k' + '/'

mitoPrefix = 'MT-'

os.system('mkdir -p ' + write_dir)

# filtering parameters
gene_number_threshold_min = 200
cell_number_threshold_min = 3
gene_number_threshold_doublet = 5000
mt_count_percent_threshold = 20
count_normalization_factor = 1e4
highly_variable_mean_min = 0.0125
highly_variable_mean_max = 3
highly_variable_dispersion_min = 0.5
stdev_threshold_max = 10
    
# clustering paramters
number_of_neighbors = 20
number_of_PCAs = 60
number_of_PCAs_neighbors = 40
random_state_variable = 1
sc.settings.verbosity = 0                                # verbosity: errors (0), warnings (1), info (2), hints (3)

#%%
toolsIndex=['CR', 'Sfu', 'Ssp', 'kb', 'Ask', 'Adf']#, 'Ase', 'Apa', 'Afu']
tools=pd.DataFrame(index=toolsIndex, columns=['name', 'mtxDir'])
tools.name['CR']='CellRanger'
tools.mtxDir['CR']=mtxPrefix + 'CellRanger_5.0.1/' +genome1+ '/standard/default/' +sample1+ '/Run1/outs/filtered_feature_bc_matrix'

tools.name['Sfu']='STAR fullSA'
tools.mtxDir['Sfu']=mtxPrefix + 'STAR_2.7.9x/' +genome1+ '/fullSA/10X_CR4_noSAM/' +sample1+ '/Solo.out/Gene/raw/'

tools.name['Ssp']='STAR sparseSA'
tools.mtxDir['Ssp']=mtxPrefix + 'STAR_2.7.9x/' +genome1+ '/sparseSA3/10X_CR4_noSAM/' +sample1+ '/Solo.out/Gene/raw/'

tools.name['kb']='kallisto|bustools'
tools.mtxDir['kb']=mtxPrefix + 'kbpy_0.25.0/' +genome1+ '/standard_1/default/' +sample1+ '/counts_unfiltered'

tools.name['Adf']='alevin-fry decoyFull'
tools.mtxDir['Adf']=mtxPrefix + 'salmon-alevin-fry_1.4.0_0.1.0/' +genome1+ '/decoyFull/rad/' +sample1+ '/gpl_knee/quant_cr-like/alevin'

tools.name['Ask']='alevin-fry sketch'
tools.mtxDir['Ask']=mtxPrefix + 'salmon-alevin-fry_1.4.0_0.1.0/' +genome1+ '/standard/sketch_rad/' +sample1+ '/gpl_knee/quant_cr-like/alevin'

for toolID in tools.index:
    print(tools.mtxDir[toolID])

#%% load mtx and save as h5. Trim barcodes to 16b (for CR)
cellsCommon = []

for toolID in tools.index:
    fOut = write_dir + toolID + '_raw.h5ad'
    print(toolID)
    if not os.path.exists(fOut):
        adata=sc.read_10x_mtx(tools.mtxDir[toolID])
        adata.obs_names = [c[0:16] for c in adata.obs_names]
        adata.write(fOut)
    else:
        adata=sc.read(write_dir + toolID + '_raw.h5ad')

    cells1 = adata.obs_names.to_list()
    print(toolID, len(cells1), len(cellsCommon))
    if len(cellsCommon)!=0:
        cellsCommon = list(set(cellsCommon) & set(cells1))
    else:
        cellsCommon = cells1
    print(toolID, len(cells1), len(cellsCommon))
    
#%% select common cells
for toolID in tools.index:
    print(toolID)
    fOut = write_dir + toolID + '_common.h5ad'
    if os.path.exists(fOut):
        continue    
    adata=sc.read(write_dir + toolID + '_raw.h5ad')
    adata=adata[cellsCommon, :]
    adata.write(fOut)

#%% filter CR and select its cell as common for all tools
for toolID in ['CR']:
    print(toolID)

    adata=sc.read(write_dir + toolID + '_common.h5ad')
    sc.pp.filter_genes(adata, min_cells=cell_number_threshold_min)
    print('MinCells',adata.X.shape)   
    sc.pp.filter_cells(adata, min_genes=gene_number_threshold_min)
    print('MinGenes',adata.X.shape)
    
    #%%
    adata.var['mt'] = adata.var_names.str.startswith(mitoPrefix)  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    sc.pl.violin(adata, ['pct_counts_mt'], jitter=0.4, xlabel=toolID, title=toolID)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', title=toolID)
    
    # remove genes with high  of MT reads
    adata = adata[adata.obs.pct_counts_mt < mt_count_percent_threshold, :]
    print('After mito removal:',adata.X.shape)
    cellsCommon1 = adata.obs_names
       
#%% commonCells1, normalize, log, HVG
for toolID in tools.index:
    print(toolID)
    fOut = write_dir + toolID + '_log.h5ad'
    if os.path.exists(fOut):
        continue
    
    adata = sc.read(write_dir + toolID + '_common.h5ad')
    
    adata = adata[cellsCommon1,:]
    print('cellsCommon1',adata.X.shape)   

    # filter genes by min_cells    
    sc.pp.filter_genes(adata, min_cells=cell_number_threshold_min)
    print('MinCells',adata.X.shape)   

    # identify highly-variable genes
    sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=3000)

    # total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, so that counts become comparable among cells
    sc.pp.normalize_total(adata, target_sum=count_normalization_factor)
    
    # logarithmize the data
    sc.pp.log1p(adata) # computes X=log(X+1), where log denotes the natural logarithm unless a different base is given
    

    adata.raw = adata
    adata.write(fOut)    
    
#%% hvg, regress, scale
for toolID in tools.index:
    print(toolID)
    fOut = write_dir + toolID + '_scale.h5ad'
    if os.path.exists(fOut):
        continue    
    
    adata=sc.read(write_dir + toolID + '_log.h5ad')
    adata=adata[:, adata.var.highly_variable]
    # regress out effects of total counts per cell and the percentage of mitochondrial genes expressed
    # scale the data to unit variance
    #adata.var['mt'] = adata.var_names.str.startswith(mitoPrefix)  # annotate the group of mitochondrial genes as 'mt'
    #sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)   
    #sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    
    # scale each gene to unit variance
    # clip values exceeding standard deviation 10
    sc.pp.scale(adata, max_value=stdev_threshold_max)
    adata.write(fOut)

#%% PCA and neighbors
for toolID in tools.index:
    print(toolID)
    fOut = write_dir + toolID + '_neighbors.h5ad'
    if os.path.exists(fOut):
        continue    
    
    adata=sc.read(write_dir + toolID + '_scale.h5ad')
    # principal component analysis
    # reduce the dimensionality of the data by running principal component analysis (PCA), 
    # which reveals the main axes of variation and denoises the data
    sc.pp.pca(adata, use_highly_variable=True, svd_solver='full', n_comps=number_of_PCAs, random_state=random_state_variable)
    
    # compute the neighborhood graph of cells using the PCA representation of the data matrix
    from neighbors_dense import neighbors_dense
    neighbors_dense(adata, n_neighbors=number_of_neighbors, n_pcs=number_of_PCAs_neighbors, use_rep='X_pca', random_state=random_state_variable, method='gauss', knn=True)
    sc.tl.umap(adata, random_state=random_state_variable)
    adata.write(fOut)
