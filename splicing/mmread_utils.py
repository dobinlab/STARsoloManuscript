from collections import Counter
import multiprocessing as mp
import numpy as np
import os
import re
from scipy.io import mmread
import scipy.sparse as sps
from scipy.stats import mannwhitneyu, ttest_ind

def adjust_mats(sparse_mats_sj, sparse_mats_gene, all_sjs, sj_annots):
    """
    Return adjusted SJ and gene mats (in CSR format), and new SJ list
    """
    new_genes_mat = []
    dup_genes_mat = []
    new_genes_idx = []
    dup_genes_idx = []
    all_gene_ids = [] # used as tie-breaker when sorting SJs at the end
    sj_dup_num = np.zeros(sparse_mats_sj.shape[0], dtype=int)
    sj_ord = np.asarray(sorted(all_sjs))
    for i,sj in enumerate(sj_ord):
        try:
            gene_ids = sj_annots[tuple(sj)]
        except KeyError:
            new_genes_mat.append(sps.coo_matrix((1,sparse_mats_gene.shape[1])))
            new_genes_idx.append(-1)
            all_gene_ids.append(sparse_mats_gene.shape[0]+1)
            continue

        all_gene_ids.extend(gene_ids)
        if len(gene_ids) > 1:
            sj_dup_num[i] = len(gene_ids) - 1

        for j,g in enumerate(sorted(gene_ids)):
            to_add = sparse_mats_gene[g-1,:]
            if j > 0:
                dup_genes_mat.append(to_add)
                dup_genes_idx.append(g-1)
            else:
                new_genes_mat.append(to_add)
                new_genes_idx.append(g-1)

    sj_dup = [np.tile(sj_ord[i], (n,1)) for i,n in enumerate(sj_dup_num)]
    sj_counts_dup = sps.vstack([sparse_mats_sj[i,:] \
        for i,n in enumerate(sj_dup_num) for _ in range(n)])
    sj_counts_dup = sps.vstack([sparse_mats_sj, sj_counts_dup], 'csr')

    sj_ord = np.concatenate((sj_ord, *sj_dup), axis=0)
    print(len(all_sjs), sj_ord.shape, len(all_gene_ids), len(new_genes_idx),
        len(new_genes_mat), len(dup_genes_idx), len(dup_genes_mat), flush=True)
    # Add gene ids for a tie-breaker for the duplicated SJs
    sort_sjs = np.concatenate((sj_ord, np.reshape(all_gene_ids, (-1,1))),
        axis=1)
    # Flip and transpose matrix so that chrom is last column
    #  (first sort col)
    sort_sjs = np.fliplr(sort_sjs).T
    sort_idx = np.lexsort(sort_sjs)

    new_sj = sj_counts_dup[sort_idx,:]
    g_idx = np.concatenate([new_genes_idx, dup_genes_idx])[sort_idx]
    new_gene = sps.vstack(new_genes_mat+dup_genes_mat, 'csr')[sort_idx,:]
    # new_gene = sps.coo_matrix((len(g_idx), sparse_mats_gene.shape[1])).tocsr()
    # nonzero_idx = g_idx != -1
    # print(new_gene.shape, sparse_mats_gene.shape, len(nonzero_idx),
    #     sum(nonzero_idx), max(g_idx), max(g_idx[nonzero_idx]), min(g_idx),
    #     min(g_idx[nonzero_idx]), g_idx[nonzero_idx].shape, flush=True)
    # indptr = sparse_mats_gene.tocsr().indptr
    # print(sum(np.diff(indptr) < 0))
    # sparse_mats_gene = sparse_mats_gene.tocsr()
    # sparse_mats_gene.indices = sparse_mats_gene.indices.astype(np.int64)
    # new_gene[nonzero_idx,:] = sparse_mats_gene[g_idx[nonzero_idx],:]
    sj_ord = sj_ord[sort_idx,:]
    return(new_sj, new_gene, sj_ord, g_idx)

def bh_test(pvals, fdr=0.05):
    """
    Use BH test to find SJs with significant difference. All SJs that failed the
    rel_filt in calc_p_vals are automatically assumed non-significant.
    """
    bad_idx = pvals == -1
    good_pvals = pvals[~bad_idx]

    sort_idx = np.argsort(good_pvals)
    pval_list = good_pvals[sort_idx]
    n_val = len(pval_list)
    bh_list = np.arange(1, n_val+1) / n_val * fdr
    try:
        # Take last index of all indices whose values are <= their BH crit value
        cutoff_idx = np.arange(1, n_val+1)[pval_list <= bh_list][-1]
    except IndexError:
        cutoff_idx = 0

    # All of the indices up to cutoff_idx are indices of significant SJs
    idx = np.zeros(len(pval_list), dtype=bool)
    idx[sort_idx[:cutoff_idx]] = True
    sig_idx = np.zeros(len(pvals), dtype=bool)
    sig_idx[~bad_idx] = idx

    return(sig_idx)

def bootstrap_iter(sj_counts, gene_counts):
    """
    Perform one iteration of the below bootstrapping. Returns the sum of sampled
    SJ counts normalized by sum of sampled gene counts.
    """
    ## Generate random index
    n_cells = len(sj_counts)
    idx = np.random.choice(range(n_cells), n_cells)

    ## Perform summing and normalization
    sj = sum(sj_counts[idx])
    g = sum(gene_counts[idx])

    if g == 0:
        u = 0
    else:
        u = sj/g

    return(u)

def bootstrap_mp(sj_counts, gene_counts, n_iter):
    """
    Multiprocessing wrapper function for performing the bootstrapping.
    """
    return([bootstrap_iter(sj_counts, gene_counts) for _ in range(n_iter)])

def bootstrap_rel_u(sj_counts, gene_counts, n_iter=100, proc=None):
    """
    Normalize SJ counts by gene counts. Returns a mat of individual cells
    normalized (in CSR), and an (n_gene,n_iter) array of the sum over
    sampled (with replacement) cells of SJ counts normalized by sum over sampled
    cells of gene counts.
    """
    ## Normalize each cell individually
    # Create a new array to hold the results
    norm_counts = sps.lil_matrix(sj_counts.shape)
    # Only do division where gene_counts are non-zero, keep zero elsewhere
    div_idx = gene_counts.nonzero()
    print(sj_counts.shape, gene_counts.shape, norm_counts.shape)
    norm_counts[div_idx] = sj_counts[div_idx] / gene_counts[div_idx]
    print(sj_counts.shape, gene_counts.shape, norm_counts.shape, flush=True)

    ## Normalize cell sums via bootstrapping
    # Create a new array to hold the results
    if proc is None:
        norm_counts_sum = np.zeros((sj_counts.shape[0],n_iter))
        for i in range(sj_counts.shape[0]):
            norm_counts_sum[i,:] = [\
                bootstrap_iter(sj_counts[i,:].toarray().flatten(),
                    gene_counts[i,:].toarray().flatten()) \
                for _ in range(n_iter)]
    else:
        with mp.Pool(proc) as pool:
            cmd_args = [(sj_counts[i,:].toarray().flatten(),
                gene_counts[i,:].toarray().flatten(), n_iter) \
                for i in range(sj_counts.shape[0])]
            norm_counts_sum = np.asarray(pool.starmap(bootstrap_mp, cmd_args))

    return(norm_counts, norm_counts_sum)

def calc_mat_n50(mat):
    """
    Calculate "N50" values for a matrix of gene counts. Exclude all genes with
    zero expression across all cells. Returns a single value.
    """
    gene_sums = np.asarray(mat.sum(axis=1)).flatten()
    gene_sums = gene_sums[gene_sums>0]

    total_expr = sum(gene_sums)
    expr_counts = Counter(gene_sums)

    expr_sum = 0
    for e in sorted(expr_counts.keys()):
        c = expr_counts[e]
        expr_sum += e * c

        if expr_sum >= (total_expr / 2):
            break

    return(e)

def calc_pair_lfc(nc1, nc2):
    """
    Calculate log2(U1/U2). Takes two norm_counts_sum results from calc_u().
    """
    ## First find places where each norm counts is non-zero
    nonzero_idx1 = nc1 != 0
    nonzero_idx2 = nc2 != 0
    try:
        nonzero_idx1 = nonzero_idx1.toarray().flatten()
        nonzero_idx2 = nonzero_idx2.toarray().flatten()
    except AttributeError:
        pass
    # Locations where both are nonzero
    both_nonzero_idx = nonzero_idx1 & nonzero_idx2
    # Locations where both are zero
    both_zero_idx = ~nonzero_idx1 & ~nonzero_idx2
    # Locations where one is zero and one is not
    diff_idx = ~both_nonzero_idx & ~both_zero_idx
    ## For each cell type, find locations where that ct is nonzero and only one
    ##  of the two cts are nonzero -> know that that ct is the nonzero one and
    ##  the other ct is zero
    nonzero_idx1 &= diff_idx
    nonzero_idx2 &= diff_idx

    ## Make results mat and use nonzero entries
    pair_lfc = np.zeros(nc1.shape)

    ## Calculate lfc for locations where both cts are nonzero
    nonzero_lfc = np.log2(nc1[both_nonzero_idx] / nc2[both_nonzero_idx])
    pair_lfc[both_nonzero_idx] = nonzero_lfc

    ## Assign values for places where at least one nc is 0
    # If only ct1 is nonzero -> +inf
    pair_lfc[nonzero_idx1] = np.inf
    # If only ct2 is nonzero -> -inf
    pair_lfc[nonzero_idx2] = -np.inf

    # pair_lfc[pair_lfc>5] = 5
    # pair_lfc[pair_lfc<-5] = -5
    
    assert sum(both_nonzero_idx) + sum(both_zero_idx) + sum(nonzero_idx1) + \
        sum(nonzero_idx2) == pair_lfc.shape[0]

    # print(np.histogram(pair_lfc.toarray().flatten(), bins=25))
    return(pair_lfc)

def calc_pair_rel_usage(nc1, nc2):
    """
    Calculate (U1-U2)/(U1+U2). Takes two norm_counts_sum results from calc_u().
    """
    ## First calculate numerator and denominator
    num = nc1 - nc2
    denom = nc1 + nc2

    ## Make results mat and use nonzero entries of denom
    pair_rel = sps.lil_matrix((denom.shape[0],1))
    div_idx = denom.nonzero()
    pair_rel[div_idx] = num[div_idx] / denom[div_idx]

    return(pair_rel.toarray().flatten())

def calc_p_vals(nc1, nc2, pair_rel, rel_filt=0.6, test='mw'):
    """
    Calculate p-values of a Mann-Whitney test for a pair of normalized counts.
    For SJs that don't pass the rel_filt on pair_rel, set p-val to -1 to keep
    the shapes the same but to denote not to include it in BH test.
    """

    idx = np.abs(pair_rel) <= rel_filt
    pvals = np.repeat(-1., len(pair_rel))
    total_0 = 0
    for i in range(len(pair_rel)):
        if idx[i]: continue
        tmp1 = nc1[i,:]
        tmp2 = nc2[i,:]
        try:
            tmp1 = tmp1.toarray().flatten()
            tmp2 = tmp2.toarray().flatten()
        except AttributeError:
            pass
        if len(set(np.concatenate((tmp1, tmp2)))) == 1:
            # print(f'SJ counts but no gene: {i}')
            continue
        if test == 'mw':
            pvals[i] = mannwhitneyu(tmp1, tmp2)[1]
        elif test == 'ttest':
            pvals[i] = ttest_ind(tmp1, tmp2, equal_var=False,
                nan_policy='omit')[1]
        # u, p = mannwhitneyu(tmp1, tmp2, alternative='two-sided')
        # if p == 0:
        #     print(f'Score and pval 0: {i}')
        #     print(u)
        #     print(len(tmp1), len(tmp2))
        #     print(len(set(tmp1)), len(set(tmp2)), flush=True)
        #     total_0 += 1
        # pvals[i] = p

        # except ValueError as e:
        #     print('1', i)
        #     print('4', idx[i])
        #     print('5', pair_rel[i])
        #     print('7', len(set(tmp1)), len(set(tmp2)),
        #         len(set(np.concatenate((tmp1, tmp2)))), flush=True)
            # raise e

    # print(f'Total 0 pvals: {sum(pvals == 0)}')
    # print(f'Counted 0: {total_0}')
    # pvals[~idx] = [mannwhitneyu(nc1[i,:].toarray()[0],
    #     nc2[i,:].toarray()[0])[1] for i in range(len(pair_rel)) if not idx[i]]
    return(pvals)

def calc_u(sj_counts, gene_counts):
    """
    Normalize SJ counts by gene counts. Returns a mat of individual cells
    normalized (in CSR), and a vector of the sum over all cells of SJ counts
    normalized by sum over all cells of gene counts.
    """
    ## Normalize each cell individually
    # Create a new array to hold the results
    norm_counts = sps.lil_matrix(sj_counts.shape)
    # Only do division where gene_counts are non-zero, keep zero elsewhere
    div_idx = gene_counts.nonzero()
    norm_counts[div_idx] = sj_counts[div_idx] / gene_counts[div_idx]

    ## Normalize cell sums
    # Create a new array to hold the results
    norm_counts_sum = sps.lil_matrix((sj_counts.shape[0],1))
    # Sum SJ counts and gene counts
    sj_sum = sj_counts.sum(axis=1)
    gene_sum = gene_counts.sum(axis=1)
    # Only do division where gene_counts are non-zero, keep zero elsewhere
    div_idx = gene_sum.nonzero()
    norm_counts_sum[div_idx] = sj_sum[div_idx] / gene_sum[div_idx]

    return(norm_counts, norm_counts_sum)

def filter_bootstrap(bootstrap, n_std=3):
    """
    Finds an index of SJs whose mean U values in bootstrap are within n_std
    standard deviations of 0 or 1.
    """
    ## Calculate means and standard devs
    b_means = np.mean(bootstrap, axis=1)
    b_stds = np.std(bootstrap, axis=1)

    ## Test values
    idx = (b_means - n_std*b_stds <= 0) | (b_means + n_std*b_stds >= 1)
    return(~idx)

def filter_gene_counts(gene_counts, gene_counts_orig, filt=25):
    """
    Filter sparse mat (in CSR format) based on gene counts for all genes.
    """
    gene_sums = np.asarray(gene_counts.sum(axis=1)).flatten()
    gene_sums_orig = np.asarray(gene_counts_orig.sum(axis=1)).flatten()

    idx = gene_sums_orig > 0
    cutoff = np.percentile(gene_sums_orig[idx], (100-filt))

    # print(sum(gene_sums_orig > cutoff))
    print(cutoff)
    return(gene_sums > cutoff)

def filter_sj_counts(sj_counts, filt=0):
    """
    Filter sparse mat (in CSR format) based on SJ counts for all SJs.
    """
    sj_sums = np.asarray(sj_counts.sum(axis=1)).flatten()

    return(sj_sums > filt)

def filter_mats(sj_counts, gene_counts, sj_filt=0, sj_ft='counts',
    g_filt=25, g_ft='perc'):
    """
    Filter sparse mats (in CSR format) based on an SJ filter and gene filter.
    The *_ft args are the filter type (either 'counts' to filter by raw counts
    or 'perc' to filter by percent).
    Returns filtered mats.
    """
    sj_sums = np.asarray(sj_counts.sum(axis=1)).flatten()
    gene_sums = np.asarray(gene_counts.sum(axis=1)).flatten()
    cell_gene_sums = np.asarray(gene_counts.sum(axis=0)).flatten()


    ## Calculate SJ and gene cutoffs (if necessary)
    if g_ft == 'perc':
        idx = gene_sums > 0
        gene_cutoff = np.percentile(gene_sums[idx], (100-g_filt))
    else:
        gene_cutoff = g_filt
    if sj_ft == 'perc':
        sj_cutoff = np.percentile(sj_sums, sj_filt)
    else:
        sj_cutoff = sj_filt

    g_idx = (sj_sums > sj_cutoff) & (gene_sums > gene_cutoff)
    # bc_idx = cell_gene_sums > 0

    return(g_idx)

def get_all_sjs(in_dir, sj_annots):
    """
    Return a set of all present SJs across all samples and a list of lists of
    tuples of SJs for each sample.

    SJ = (chrom, start_pos, end_pos, strand)
    """
    samp_dirs = [d for d in os.listdir(in_dir) \
        if os.path.isdir(f'{in_dir}/{d}')]

    fns = [f'{in_dir}/{d}/SJ.out.tab' for d in samp_dirs]

    sj_lists = []
    all_sjs = set()
    for fn in fns:
        sjs = []
        for sj in open(fn, 'r'):
            sj = sj.strip().split('\t')

            if sj[5] == '0':
                # unannotated
                continue

            sj[0] = sj[0].strip('chr')
            if sj[0] == 'X':
                sj[0] = '23'
            elif sj[0] == 'Y':
                sj[0] = '24'
            elif sj[0][0] in {'M', 'm'}:
                sj[0] = '25'
            elif not re.match(r'^[0-9]+$', sj[0]):
                continue

            signs = []
            # if sj[3] == '1':
            #     signs.append(1)
            # elif sj[3] == '2':
            #     signs.append(2)
            if sj[3] == '1' or sj[3] == '2':
                signs.append(sj[3])
            else:
                for s in ['1', '2']:
                    if tuple(sj[:3] + [s]) in sj_annots:
                        signs.append(s)
                # if tuple(sj[:3] + ['1']) in sj_annots:
                #     signs.append('1')
                # if tuple(sj[:3] + ['2']) in sj_annots:
                #     signs.append('2')

                if len(signs) == 0:
                    signs.append(-1)

            new_sjs = [tuple(list(map(int, sj[:3])) + [s]) for s in signs]
            sjs.append(tuple(new_sjs))
            all_sjs.update(new_sjs)

        sj_lists.append(sjs)

    return(all_sjs, sj_lists)


def load_all_gene_counts(in_dir, counts_type='raw'):
    """
    Load all matrix gene count files and combine them into one large sparse
    array in CSC format.

    Returns the large CSC matrix and the list of matrix file names.
    """
    samp_dirs = [d for d in os.listdir(in_dir) \
        if os.path.isdir(f'{in_dir}/{d}')]

    fns = [f'{in_dir}/{d}/Solo.out/Gene/{counts_type}/matrix.mtx' \
        for d in samp_dirs]

    # mats = []
    # for fn in fns:
    #     print(fn, flush=True)
    #     m = mmread(fn)
    #     print(m.shape, flush=True)
    #     mats.append(m)
    return(sps.hstack([mmread(fn) for fn in fns], 'csc', int), fns)

def load_all_sj_counts(in_dir, all_sjs, sj_lists):
    """
    Load all matrix gene count files and combine them into one large sparse
    array in CSC format. Need to make sure that all mats have the same SJs in
    the same order so we can combine them together.

    Returns the large CSC matrix and the list of matrix file names.
    """
    """
    (will want mats in CSR format for this)
    For each sample, make sure that all of the all_sjs are present.
    Then, remove any extra rows
    Finally, sort rows by np.argsort(all_sjs) (make sure to use the order that
                                                the SJs are actually in for that
                                                sample)
    """

    samp_dirs = [d for d in os.listdir(in_dir) \
        if os.path.isdir(f'{in_dir}/{d}')]
    fns = [f'{in_dir}/{d}/Solo.out/SJ/raw/matrix.mtx' \
        for d in samp_dirs]
    # all_counts = sps.hstack([mmread(fn) for fn in fns], 'csr', int)

    cmd_args = [(fn, sj_lists[i], all_sjs) for i,fn in enumerate(fns)]
    n_procs = min(30, len(fns))
    with mp.Pool(n_procs) as pool:
        all_counts = pool.starmap(load_sj_mp, cmd_args)

    return(sps.hstack(all_counts, 'csc', int), fns)

def load_sj_mp(fn, sj_list, all_sjs):
    """
    Function to use for multiprocessing
    """
    ## Find which of the rows already present we'll keep
    idx = np.asarray([np.any([sj in all_sjs for sj in sjs]) \
        for sjs in sj_list])
    ## Find all SJs present in this mat
    cur_sjs = set(sj_list[0])
    for sjs in sj_list[1:]:
        cur_sjs.update(sjs)
    ## Find which SJs we need to add
    to_add = all_sjs - cur_sjs

    ## Find duplicate SJs to add
    ## Duplicates can exist for cases where the sign isn't listed for the
    ##  discovered SJ, and SJs exist in genome on both strands
    dup_idx = np.asarray([np.sum([sj in all_sjs for sj in sjs]) > 1 \
        for sjs in sj_list])
    print(f'{fn}: {sum(dup_idx)}', flush=True)
    dup_to_add = [sjs[1] for j,sjs in enumerate(sj_list) if dup_idx[j]]

    assert len(to_add) + len(dup_to_add) + sum(idx) == len(all_sjs)

    ## Add all new bools
    idx = np.concatenate([idx, [True]*(len(to_add)+len(dup_to_add))])

    ## Load matrix
    m = mmread(fn).tocsr()
    ## Get dup vals
    dup_mat = m[dup_idx,:]
    ## Make matrix of zeros to add at the bottom
    to_add_mat = sps.coo_matrix((len(to_add), m.shape[1]), int)
    ## Combine mats
    m = sps.vstack([m, to_add_mat, dup_mat], 'csr', int)

    ## Take approprate subset of the array
    m = m[idx,:]

    ## Sort SJs (use first SJ for each row that's in all_sjs)
    sort_sjs = [sjs[0] if sjs[0] in all_sjs else sjs[1] \
        for j,sjs in enumerate(sj_list) if idx[j]]
    # Doesn't matter what order to_add is in since they're all zeros
    # Ultimately end up sorting all_sjs
    sort_sjs = np.asarray(sort_sjs + list(to_add) + dup_to_add)
    # Flip and transpose matrix so that chrom is last column
    #  (first sort col)
    sort_sjs = np.fliplr(sort_sjs).T
    sort_idx = np.lexsort(sort_sjs)
    print(f'{fn}: {sort_idx.shape}', flush=True)
    m = m[sort_idx,:]

    return(m)

def split_counts_mat(all_counts, mtx_fns, annots_dict, cts=None):
    """
    Create dictionary of cell type -> CSR matrix (combined across all samples)
    """
    if cts is None:
        cts = {v for v in annots_dict.values()}

    cts = sorted(cts)

    bcs = []
    for fn in mtx_fns:
        bc_fn = f'{os.path.dirname(fn)}/barcodes.tsv'
        sample = re.split(r'/+', fn)[-5]
        bcs.extend([f'{sample}_{line.strip()}' for line in open(bc_fn, 'r')])
    bcs = np.asarray(bcs)

    cell_cts = np.asarray([annots_dict[bc] if bc in annots_dict else np.nan \
        for bc in bcs])
    sparse_mats = {}
    cell_bcs = {}
    for ct in cts:
        # idx = np.asarray([bc in annots_dict and annots_dict[bc] == ct \
        #     for bc in bcs])
        idx = cell_cts == ct
        sparse_mats[ct] = all_counts.tocsc()[:,idx].tocsr()
        cell_bcs[ct] = bcs[idx]

    return(sparse_mats, cell_bcs)