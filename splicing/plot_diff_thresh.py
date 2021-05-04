import argparse
import itertools as it
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import mmread_utils as mmu
import numpy as np
import pandas
import pickle as pkl
import re
import seaborn as sns

def load_stats_dict(fn):
    df = pandas.read_csv(fn, index_col=0)
    stats_dict = {}
    cts = np.unique(df['ct'])

    for ct in cts:
        idx = df['ct'].to_numpy() == ct
        stats_dict[ct] = df.iloc[idx,:-1].to_numpy()

    return(stats_dict)

def plot_all_agg(tiss_pair_stats, tiss_sig_idxs, all_thresh, fn_out, cutoff=5):
    fig, ax = plt.subplots(figsize=(18,9))

    bins = np.linspace(0, cutoff, 50)
    colors = sns.color_palette()[:3]
    total_sjs = 0
    for i,sj_filt in enumerate(all_thresh):
        idx = np.concatenate([tiss[sj_filt] for tiss in tiss_sig_idxs.values()])

        stat_max = np.concatenate([np.nanmax(np.abs(stat[sj_filt]), axis=1) \
            for stat in tiss_pair_stats.values()])
        nan_idx = np.isnan(stat_max)

        stat_max = stat_max[idx & ~nan_idx]

        stat_max[stat_max>cutoff] = cutoff
        sns.histplot(x=stat_max, bins=bins, ax=ax, legend=False, zorder=i+25,
            common_norm=False, fill=True, alpha=1, color=colors[i],
            label=sj_filt)

        total_sjs += len(stat_max)

    ax.legend(title='SJ Counts Filter')
    ax.set_xlabel(r'$\left|\mathrm{Log_2\ Fold\ Change}\right|$')
    ax.get_yaxis().set_major_formatter(
        FuncFormatter(lambda x, p: f'{int(x):,}'))

    ax.grid(True, axis='y')

    print(f'Total SJs agg: {total_sjs}', flush=True)

    fig.savefig(fn_out, bbox_inches='tight', dpi=200)
    plt.close(fig)

def plot_agg(total_pair_stats, total_sig_idxs, fn_out, cutoff=5):
    fig, ax = plt.subplots(figsize=(18,9))

    bins = np.linspace(0, cutoff, 50)
    colors = sns.color_palette()[:3]
    total_sjs = 0
    for i,(sj_filt,agg_stats) in enumerate(total_pair_stats.items()):
        idx = total_sig_idxs[sj_filt]

        stat_max = np.nanmax(np.abs(agg_stats), axis=1)
        nan_idx = np.isnan(stat_max)

        stat_max = stat_max[idx & ~nan_idx]

        stat_max[stat_max>cutoff] = cutoff
        sns.histplot(x=stat_max, bins=bins, ax=ax, legend=False, zorder=i+25,
            common_norm=False, fill=True, alpha=1, color=colors[i],
            label=sj_filt)

        total_sjs += len(stat_max)

    ax.legend(title='SJ Counts Filter')
    ax.set_xlabel(r'$\left|\mathrm{Log_2\ Fold\ Change}\right|$')
    ax.get_yaxis().set_major_formatter(
        FuncFormatter(lambda x, p: f'{int(x):,}'))

    ax.grid(True, axis='y')

    print(f'Total SJs agg: {total_sjs}', flush=True)

    fig.savefig(fn_out, bbox_inches='tight', dpi=200)
    plt.close(fig)

def plot_hists(stat_df, fn_out, cutoff=5):
    fig, ax = plt.subplots(figsize=(18,12))

    ## Set values that exceed the cutoff as the cutoff value
    stat_df.loc[stat_df['stat_vals'] > cutoff, 'stat_vals'] = cutoff
    stat_df.loc[stat_df['stat_vals'] < -cutoff, 'stat_vals'] = -cutoff
    # stat_df.loc[:,'thresh'] = stat_df['thresh'].astype(str)

    bins = np.linspace(-cutoff, cutoff, 50)
    threshes = sorted(np.unique(stat_df['thresh']))
    colors = sns.color_palette()[:3]
    print(colors)

    total_sjs = 0
    for i,t in enumerate(threshes):
        sns.histplot(data=stat_df.loc[stat_df['thresh']==t,:], x='stat_vals',
            bins=bins, multiple='layer', ax=ax, legend=False, zorder=i+25,
            common_norm=False, fill=True, alpha=1, color=colors[i], label=t)

        total_sjs += sum(stat_df['thresh']==t)
    # sns.histplot(data=stat_df, x='stat_vals', hue='thresh', bins=bins,
    #     multiple='layer', ax=ax, legend='full', common_norm=False, fill=True,
    #     palette='tab10', alpha=1)
    # element='step', cumulative=True, stat='density', 
    # sns.kdeplot(data=stat_df, x='stat_vals', hue='thresh', cumulative=True,
    #     multiple='layer', ax=ax)

    ax.legend()
    ax.set_xlabel(r'$\mathrm{Log_2\ Fold\ Change}$')
    ax.get_yaxis().set_major_formatter(
        FuncFormatter(lambda x, p: f'{int(x):,}'))
    # handles, labels = ax.get_legend_handles_labels()
    # ax.legend(handles=handles, labels=labels, title='SJ Count Threshold')

    ax.grid(True, axis='y')

    print(f'Total SJs {fn_out}: {total_sjs}')

    fig.savefig(fn_out, bbox_inches='tight', dpi=200)
    plt.close(fig)

def plot_scatter(stat_df, stats_dict, ct1, ct2, fn_out, sj_ord, sj_dict,
    genes_dict):
    # fig, axes = plt.subplots(nrows=3, figsize=(18,36))
    fig, ax = plt.subplots(figsize=(18,12))

    u1 = np.median(stats_dict[ct1], axis=1).reshape([-1,1])
    u2 = np.median(stats_dict[ct2], axis=1).reshape([-1,1])

    sj_df = None
    for i,thresh in enumerate(np.unique(stat_df['thresh'])):
        print(i, thresh, flush=True)
        # ax = axes.flatten()[i]

        temp_u = np.concatenate([u1, u2], axis=1)
        t_idx = stat_df['thresh'].to_numpy()==thresh
        idx = stat_df.loc[t_idx,'idx'].to_numpy()
        temp_u = temp_u[idx]

        u_max = np.max(temp_u, axis=1)
        lfc = stat_df.loc[t_idx,'stat_vals'].to_numpy()[idx]

        pos_inf_idx = lfc == np.inf
        neg_inf_idx = lfc == -np.inf
        inf_idx = pos_inf_idx | neg_inf_idx

        # sns.scatterplot(x=np.max(temp_u, axis=1),
        #     y=stat_df.loc[t_idx,'stat_vals'].to_numpy()[idx], ax=ax)
        sns.scatterplot(x=u_max[~inf_idx], y=lfc[~inf_idx],
            color=sns.color_palette()[i], label=thresh, ax=ax)
        sns.scatterplot(x=u_max[pos_inf_idx],
            y=[max(lfc[~inf_idx])]*sum(pos_inf_idx),
            color=sns.color_palette()[i], label=thresh, marker='^', ax=ax)
        sns.scatterplot(x=u_max[neg_inf_idx],
            y=[min(lfc[~inf_idx])]*sum(neg_inf_idx),
            color=sns.color_palette()[i], label=thresh, marker='v', ax=ax)

        # sig_idx = (np.max(temp_u, axis=1) > 0.6) & \
        #     (np.abs(stat_df.loc[t_idx,'stat_vals'].to_numpy()[idx]) > 2)
        # print('^^', ct1, ct2, thresh, sum(sig_idx))
        # for sj in sj_ord[idx,:][sig_idx,:]:
        #     print('^^', sj, flush=True)

        df = pandas.DataFrame(sj_ord[idx,:])
        df['thresh'] = thresh
        # df['u_max'] = u_max
        df['u1'] = temp_u[:,0]
        df['u2'] = temp_u[:,1]
        df['lfc'] = lfc

        ## Add SJ count info
        sj1 = np.asarray(sj_dict[ct1].sum(axis=1)).flatten()
        sj2 = np.asarray(sj_dict[ct2].sum(axis=1)).flatten()

        sj1 = sj1[idx]
        sj2 = sj2[idx]

        df['sj1'] = sj1
        df['sj2'] = sj2

        gene1 = np.asarray(genes_dict[ct1].sum(axis=1)).flatten()
        gene2 = np.asarray(genes_dict[ct2].sum(axis=1)).flatten()

        gene1 = gene1[idx]
        gene2 = gene2[idx]

        df['gene1'] = gene1
        df['gene2'] = gene2

        if sj_df is None:
            sj_df = df.copy()
        else:
            sj_df = pandas.concat([sj_df, df.copy()])

        # ax.set_xlabel('Max U')
        # ax.set_ylabel(f'{ct1}/{ct2}')
        # ax.set_title(f'SJ Threshold {thresh}')

    ax.set_xlabel('Max U')
    # ax.set_ylabel(f'lfc({ct1}/{ct2})')
    ct1_trim = re.sub('_', r'\_', ct1.split(':')[1])
    ct2_trim = re.sub('_', r'\_', ct2.split(':')[1])
    ax.set_ylabel((r'$log_2\left(\frac{\mathrm{' f'{ct1_trim}'
        r'}}{\mathrm{' f'{ct2_trim}' r'}}\right)$'))
    h, l = ax.get_legend_handles_labels()
    h = np.asarray(h)
    l = np.asarray(l)

    ## Get indices to use for labels and handles
    unique_vals = sorted(np.unique(l), key=int)
    h = [h[l==thresh][0] for thresh in unique_vals]
    l = [l[l==thresh][0] for thresh in unique_vals]

    # if len(l) > 3:
    #     h = [h[0], h[3], h[5]]
    #     l = [l[0], l[3], l[5]]
    try:
        ax.legend(handles=h, labels=l, title='SJ Threshold', markerscale=3)
    except IndexError as e:
        print(l, h)
        raise e

    ax.grid(True)

    # ax.legend(handles=h[:3], labels=l[:3], title='SJ Threshold')
    # ax.set_title(f'SJ Threshold {thresh}')
    fig.savefig(fn_out, bbox_inches='tight', dpi=200)
    plt.close(fig)

    return(sj_df)

def run_stats(stats_dict, sj_dict, ct1, ct2, sj_filt, gene_filters, stat='lfc'):
    nc1 = stats_dict[ct1]
    nc2 = stats_dict[ct2]

    nc_sum1 = np.median(nc1, axis=1)
    nc_sum2 = np.median(nc2, axis=1)

    ## Filter based on gene counts
    # Make sure both have positive gene counts
    idx = gene_filters[ct1] & gene_filters[ct2]
    idx &= (mmu.filter_sj_counts(sj_dict[ct1], sj_filt) | \
        mmu.filter_sj_counts(sj_dict[ct2], sj_filt))
    print(f'SJs passing filter ({ct1}, {ct2}): {sum(idx)}', flush=True)

    print('Calculating stats...', flush=True)
    if stat == 'rel':
        stat_fn = mmu.calc_pair_rel_usage
    elif stat == 'lfc':
        stat_fn = mmu.calc_pair_lfc

    print(f'Using {stat_fn}', flush=True)
    stat_vals = stat_fn(nc_sum1, nc_sum2)

    print(f'inf fraction: {sum(np.abs(stat_vals[idx]) == np.inf)}/{sum(idx)}')
    stat_vals[~idx] = np.nan

    return(stat_vals, idx)

################################################################################
def get_args():
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-sd')
    parser.add_argument('-pkl')
    parser.add_argument('-t', nargs='+')
    parser.add_argument('-o')
    parser.add_argument('-df_o')

    parser.add_argument('-ct1')
    parser.add_argument('-ct2')
    parser.add_argument('-all', action='store_true')
    parser.add_argument('-agg', action='store_true')
    parser.add_argument('-agg_tiss', action='store_true')

    parser.add_argument('-g_filt', default=25, type=int)

    parser.add_argument('-pt', default='hist')

    return(parser.parse_args())

def main():
    args = get_args()

    stats_dict = load_stats_dict(args.sd)
    sj_dict, genes_dict, sj_ord, genes_dict_orig, _, _ = pkl.load(
        open(args.pkl, 'rb'))

    mpl.rcParams.update({'font.size': 30})

    if args.all:
        ## Remove cell types that have < 200 cells
        print(f'num clusters before trim: {len(sj_dict)}')
        use_cts = {ct for ct,m in sj_dict.items() if m.shape[1] > 200}
        sj_dict = {k: v for k,v in sj_dict.items() if k in use_cts}
        genes_dict = {k: v for k,v in genes_dict.items() if k in use_cts}
        genes_dict_orig = {k: v for k,v in genes_dict_orig.items() \
            if k in use_cts}
        print(f'num clusters after trim: {len(sj_dict)}', flush=True)
        
        gene_filters = {ct: mmu.filter_gene_counts(genes_dict[ct],
            genes_dict_orig[ct], args.g_filt) for ct in use_cts}

        tiss_dict = it.groupby(sorted(sj_dict.keys()),
            key=lambda x: x.split(':')[0])
        tiss_dict = {k: list(v) for k,v in tiss_dict}
        print(tiss_dict)
        ct_pairs = [p for ct_list in tiss_dict.values() \
            for p in it.combinations(ct_list, 2)]
    else:
        gene_filters = {ct: mmu.filter_gene_counts(genes_dict[ct],
            genes_dict_orig[ct], args.g_filt) for ct in [args.ct1, args.ct2]}

        ct_pairs = [(args.ct1, args.ct2)]

    ct_dict = {ct: len([p for p in it.combinations(ct_list, 2)]) \
        for ct,ct_list in tiss_dict.items()}
    print(f'{len(ct_pairs)} total comparisons', flush=True)
    print(ct_dict, flush=True)
    return

    if args.agg or args.agg_tiss:
        tiss_pair_stats = {}
        tiss_sig_idxs = {}

    for (ct1,ct2) in ct_pairs:
        stat_df = None
        for sj_filt in args.t:
            sj_filt = int(sj_filt)
            stat_vals, idx = run_stats(stats_dict, sj_dict, ct1, ct2,
                sj_filt, gene_filters)

            df = pandas.DataFrame({'stat_vals': stat_vals,
                'thresh': sj_filt, 'idx': idx})

            if stat_df is None:
                stat_df = df.copy()
            else:
                stat_df = pandas.concat([stat_df, df.copy()])

            if args.agg or args.agg_tiss:
                stat_vals = stat_vals.reshape((-1,1))
                tiss = ct1.split(':')[0]
                try:
                    stats_d = tiss_pair_stats[tiss]
                    sig_idx_d = tiss_sig_idxs[tiss]
                except KeyError:
                    stats_d = tiss_pair_stats[tiss] = {}
                    sig_idx_d = tiss_sig_idxs[tiss] = {}

                try:
                    stats_d[sj_filt] = np.concatenate(
                        (stats_d[sj_filt], stat_vals), axis=1)
                    sig_idx_d[sj_filt] = sig_idx_d[sj_filt] | idx
                except KeyError:
                    stats_d[sj_filt] = stat_vals
                    sig_idx_d[sj_filt] = idx

        fn_out = args.o.format(ct1, ct2)
        if args.pt == 'hist':
            plot_hists(stat_df, fn_out)
        elif args.pt == 'scatter':
            sj_df = plot_scatter(stat_df, stats_dict, ct1, ct2, fn_out, sj_ord,
                sj_dict, genes_dict)

            if args.df_o:
                sj_df.to_csv(args.df_o.format(ct1, ct2), na_rep='NaN')

    if args.agg:
        for tiss, t_dict in tiss_pair_stats.items():
            for sj_filt in t_dict.keys():
                print(tiss, sj_filt)
        fn_out = args.o.format('agg', 'all')
        plot_all_agg(tiss_pair_stats, tiss_sig_idxs, list(map(int, args.t)),
            fn_out)

    if args.agg_tiss:
        for tiss in tiss_pair_stats.keys():
            stats_d = tiss_pair_stats[tiss]
            sig_idx_d = tiss_sig_idxs[tiss]

            print(f'Plotting agg for tissue {tiss}', flush=True)
            fn_out = args.o.format('agg', tiss)
            plot_agg(stats_d, sig_idx_d, fn_out)


            


if __name__ == '__main__':
    main()