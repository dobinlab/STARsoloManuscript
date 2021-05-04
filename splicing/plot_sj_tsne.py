import argparse
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas
import pickle as pkl
import scipy.sparse as sps
import seaborn as sns

def get_sj_idx(sj_fn, sj_ord):
    sj_ord_idx_dict = {tuple(sj): i for i,sj in enumerate(sj_ord)}
    sj_idx = []
    for i,line in enumerate(open(sj_fn, 'r')):
        if i == 0:
            continue

        line = line.strip().split(',')
        sj = tuple(map(int, line[1:5]))
        try:
            sj_idx.append(sj_ord_idx_dict[sj])
        except KeyError:
            pass

    return(np.asarray(sj_idx))

def get_tsne_coords(annots_fn):
    annots_df = pandas.read_csv(annots_fn)
    tsne_dict = {c['cell']: (c['tissue_tSNE_1'], c['tissue_tSNE_2']) \
        for _,c in annots_df.iterrows()}

    return(tsne_dict)

def get_unique_idx(sj_ord):
    dup_idx = []
    seen_sjs = set()
    for sj in sj_ord:
        sj = tuple(sj)
        dup_idx.append(sj in seen_sjs)
        seen_sjs.add(sj)

    return(~np.asarray(dup_idx))

def make_tsne_df(tsne_dict, all_bcs):
    coords = [tsne_dict[bc] for bc_list in all_bcs.values() for bc in bc_list]
    cell_types = [ct for ct, bc_list in all_bcs.items() \
        for _ in range(len(bc_list))]
    tsne_df = pandas.DataFrame(coords, columns=['tSNE_1', 'tSNE_2'])
    tsne_df['cell_type'] = cell_types

    return(tsne_df)

def plot_tsne(tsne_df, all_u, out_base):
    mpl.rcParams.update({'font.size': 30})

    ## Adjust the color map to have gray at the beginning (224, 224, 224, 1)
    gray_val = 224 / 255
    cmap = cm.get_cmap('OrRd')
    cmap._segmentdata['red'][0] = (0., gray_val, gray_val)
    cmap._segmentdata['green'][0] = (0., gray_val, gray_val)
    cmap._segmentdata['blue'][0] = (0., gray_val, gray_val)

    for i in range(all_u.shape[0]):
        fig, ax = plt.subplots(figsize=(12,12))
        sj_exp = np.asarray(all_u[i,:]).flatten()
        sns.scatterplot(x='tSNE_1', y='tSNE_2', data=tsne_df, hue=sj_exp,
            palette=cmap, hue_norm=(0,1), s=35, ax=ax)

        ## Use AxesDivider to separate color bar from the actual plot
        div = make_axes_locatable(ax)
        cbar_ax = div.append_axes('right', size='5%', pad=0.1)

        ## Make the colorbar to replace the bad legend that seaborn makes
        norm = plt.Normalize(0, 1)
        # norm = plt.Normalize(min(sj_exp), max(sj_exp))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])

        ## Replace legend with colorbar
        ax.get_legend().remove()
        plt.colorbar(sm, cax=cbar_ax)

        fig.savefig(out_base.format(i), dpi=200, bbox_inches='tight')
        plt.close(fig)

    ## Plot the cell types
    fig, ax = plt.subplots(figsize=(12,12))
    sns.scatterplot(x='tSNE_1', y='tSNE_2', data=tsne_df, hue='cell_type',
        s=35, ax=ax)
    fig.savefig(out_base.format('ct'), dpi=200, bbox_inches='tight')
    plt.close(fig)


################################################################################
def get_args():
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-pkl')
    parser.add_argument('-sj')
    parser.add_argument('-annot')
    parser.add_argument('-o')

    parser.add_argument('-cts', nargs='+')

    return(parser.parse_args())

def main():
    args = get_args()

    sj_dict, genes_dict, sj_ord, genes_dict_orig, _, cell_bcs= pkl.load(
        open(args.pkl, 'rb'))
    
    ## Get tSNE coordinates from annots file
    tsne_dict = get_tsne_coords(args.annot)

    ## Subset sj_ord to only have unique entries
    unique_sj_idx = get_unique_idx(sj_ord)
    sj_ord = sj_ord[unique_sj_idx,:]

    ## Find the indices in the sj mat for the SJs we're looking at
    sj_idx = get_sj_idx(args.sj, sj_ord)
    print(len(set(sj_idx)))
    print(len(sj_idx), 'SJs', flush=True)

    ## Combine data across all cell types we're using
    all_bcs = {ct: cell_bcs[ct] for ct in args.cts}
    all_sj = sps.hstack([sj_dict[ct] for ct in args.cts], 'csr')
    all_genes = sps.hstack([genes_dict[ct] for ct in args.cts], 'csr')

    ## Calculate and subset U mat
    all_u = all_sj / all_genes
    all_u = all_u[unique_sj_idx,:][sj_idx,:]

    ## Plot each SJ in its own file
    tsne_df = make_tsne_df(tsne_dict, all_bcs)
    plot_tsne(tsne_df, all_u, args.o)

    ## Save the SJs
    sjs = pandas.DataFrame(sj_ord[sj_idx,:])
    sjs.to_csv(args.o.format('sjs'))

if __name__ == '__main__':
    main()