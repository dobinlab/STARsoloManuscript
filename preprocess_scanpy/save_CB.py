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

#%%
mtxPrefix='/scratch/dobin/STAR/STARsoloPreprint/count/'

toolsIndex=['CR', 'Sfu', 'Ssp', 'kb', 'Ask', 'Adf']#, 'Ase', 'Apa', 'Afu']
tools=pd.DataFrame(index=toolsIndex, columns=['name', 'mtxDir'])
tools.name['CR']='CellRanger'
tools.mtxDir['CR']=mtxPrefix + 'CellRanger_5.0.1/human_CR_3.0.0/standard/default/10X/3/pbmc_5k/20/b02/Run1/outs/filtered_feature_bc_matrix'

tools.name['Sfu']='STAR fullSA'
tools.mtxDir['Sfu']=mtxPrefix + 'STAR_2.7.8a/human_CR_3.0.0/fullSA/10X_CR4_noSAM/10X/3/pbmc_5k/20/b02/Solo.out/Gene/raw/'

tools.name['Ssp']='STAR sparseSA'
tools.mtxDir['Ssp']=mtxPrefix + 'STAR_2.7.8a/human_CR_3.0.0/sparseSA3/10X_CR4_noSAM/10X/3/pbmc_5k/20/b02/Solo.out/Gene/raw/'

tools.name['kb']='kallisto|bustools'
tools.mtxDir['kb']=mtxPrefix + 'kbpy_0.25.0/human_CR_3.0.0/standard_1/default/10X/3/pbmc_5k/20/b02/counts_unfiltered'

tools.name['Adf']='alevin-fry decoyFull'
tools.mtxDir['Adf']=mtxPrefix + 'salmon-alevin-fry_1.4.0_0.1.0/human_CR_3.0.0/decoyFull/rad_knee/10X/3/pbmc_5k/20/b02/quant_cr-like/alevin'

tools.name['Ask']='alevin-fry sketch'
tools.mtxDir['Ask']=mtxPrefix + 'salmon-alevin-fry_1.4.0_0.1.0/human_CR_3.0.0/standard/sketch_rad_knee/10X/3/pbmc_5k/20/b02/quant_cr-like/alevin'



#%% write cell barcodes for each tool
for dir_out in ['./write_pbmc5k/']:
    for toolID in tools.index:
        adata=sc.read(dir_out + toolID + '_neighbors.h5ad')
        np.savetxt(dir_out + toolID + '_CB.txt', np.array(adata.obs_names.to_list()), fmt='%s')
        np.savetxt(dir_out + toolID + '_hvg.txt', np.array(adata.var_names.to_list()), fmt='%s')
