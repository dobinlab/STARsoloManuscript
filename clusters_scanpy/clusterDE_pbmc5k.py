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

mpl.rcParams['figure.dpi'] = 150

#%%

nCl=9
sample1 = 'pbmc5k'
dirPre = '../preprocess_scanpy/write_' + sample1 + '/'
dirLeiden = '../leiden_scanpy/write_' + sample1 +'/'

dirOut1='./write_' + sample1 +'/'
dirOut =  './write_' + sample1 + '/CRclusters9/'

dirFig = './figures_' + sample1 + '/CRclusters9/'

os.system('mkdir -p ' + dirOut)
os.system('mkdir -p ' + dirFig)
#%%
sc.settings.verbosity = 3                                # verbosity: errors (0), warnings (1), info (2), hints (3)

#scanpy.set_figure_params(scanpy=True, dpi=80, dpi_save=150, frameon=True, vector_friendly=True, fontsize=14, figsize=None, color_map=None, format='pdf', facecolor=None, transparent=False, ipython_format='png2x')
sc.set_figure_params(dpi=150, dpi_save=300, fontsize=20, figsize=[8,8])
#%%
mtxPrefix='/scratch/dobin/STAR/STARsoloPreprint/maia1/count/'

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

#%%
marker_seurat = {
    'Naive CD4+ T': ['IL7R', 'CCR7'],
    'Memory CD4+':  ['IL7R', 'S100A4'],    
    'CD14+ Mono':   ['CD14', 'LYZ'],
    'B':            ['MS4A1'],
    'NK':           ['GNLY', 'NKG7'],    
    'CD8+ T':       ['CD8A'],
    'FCGR3A+ Mono': ['FCGR3A', 'MS4A7'],    
    'DC':           ['FCER1A', 'CST3'],
    'Platelet':     ['PPBP']
    }

#%% load CR
acr = sc.read(dirPre + 'CR_common.h5ad')
np.savetxt(dirOut + 'genesListFull.txt', np.asarray(acr.var_names.to_list(), dtype=str), fmt='%s')

#%% cluster CR and annotate clusters
clkey = 'CRleiden9'
toolID='CR'

acr = sc.read(dirPre + 'CR_neighbors.h5ad')
sc.tl.leiden(acr, resolution=0.36, random_state=1, key_added=clkey)

marker_seurat_types=list(marker_seurat.keys())
cluster_types = marker_seurat_types.copy()
np.savetxt(dirOut + 'clusterTypes.txt', np.asarray(cluster_types, dtype=str), fmt='%s')
acr.rename_categories(clkey, cluster_types)

#%%
plt.rc_context({'figure.figsize': (8, 8)})
sc.pl.umap(acr,color=clkey, title=tools.name[toolID] + ': Leiden clusters', save= clkey + '_' + toolID + '.pdf', legend_loc='upper right')

# on data, does not look good
#sc.pl.umap(acr,color=clkey, title=tools.name[toolID] + 'Leiden clusters', save= clkey + '_' + toolID + '_labelsOn.pdf', legend_loc='on data' )
sc.pl.stacked_violin(acr, marker_seurat, use_raw=True, groupby=clkey, figsize=[10,8], save= clkey + '_' + toolID + '.pdf' )

#%% rgg
for toolID in tools.index:
    print(toolID)
    fOut = dirOut + toolID + '_wilcoxon.h5ad'
    if os.path.exists(fOut):
        continue  
    
    adata = sc.read(dirPre + toolID + '_neighbors.h5ad')
    adata.obs[clkey] = acr.obs[clkey][adata.obs_names] # assign labels from CR
    
    print(np.count_nonzero(adata.obs_names!=acr.obs_names))
    
    for rgg_method in ['wilcoxon']: #, 'logreg']: #logreg did not work out of the box: lbfgs failed to converge (status=1): STOP: TOTAL NO. of ITERATIONS REACHED LIMIT.
        rgg_key = 'rgg_' +rgg_method
        sc.tl.rank_genes_groups(adata, clkey, method=rgg_method, tie_correct=True, key_added=rgg_key)
        #
        rgg1 = adata.uns[rgg_key]
        groups = rgg1['names'].dtype.names
        table1=pd.DataFrame(
            {group + '_' + key: rgg1[key][group]
             for group in groups for key in ['names', 'scores', 'pvals_adj', 'logfoldchanges']})
        
        print(table1)
        table1.to_csv(dirOut + toolID + '.' + rgg_key + '.csv', index=False)
        
    adata.write(fOut)
    
#%% compare SNN graph
for toolID in tools.index:
    print(toolID)
    adata = sc.read(dirPre + toolID + '_neighbors.h5ad')
    conn1 = adata.obsp['connectivities']
    conn1.maxprint = conn1.nnz
    with open(dirOut1 + toolID + '_conn.txt','w') as conn_file:
        for i in range(conn1.shape[0]):
            for j in conn1[i].nonzero()[1]:
                conn_file.write(str(i)+' ' +str(j)+' '+str(conn1[i,j])+'\n')
        conn_file.close() 