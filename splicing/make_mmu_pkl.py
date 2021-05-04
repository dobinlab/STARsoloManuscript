import argparse
import mmread_utils as mmu
import pickle as pkl

def parse_annots(fn):
    """
    Columns are:
    * 0 - run_barcode
    * 1 - cell_ontology_class
    * 17 - tissue
    """

    annots_dict = {}
    for i, line in enumerate(open(fn)):
        if i == 0: continue

        line = line.split(',')
        if len(line[1]) == 0:
            continue
        ## Combine cell types with spaces in them
        elif line[1][0] == '"':
            annot = [line[1]]
            idx = 2
            while line[idx][-1] != '"':
                annot.append(line[idx])
                idx += 1
            annot.append(line[idx])
            annot = ';'.join(annot)
            annot = annot[1:-1]
        else:
            annot = line[1]
            idx = 1

        re_rep = r"\s+"
        cell_type = (f'{re.subn(re_rep, "_", line[idx+16])[0]}:'
        f'{re.subn(re_rep, "_", annot)[0]}')
        annots_dict[line[0]] = cell_type

    return(annots_dict)

def load_sj_genes_num(gd):
    """
    Return dictionary structure:
    sj (chrom, start_pos, end_pos, strand):
        {gene_ids}
    """
    sj_list = f'{gd}/sjdbList.fromGTF.out.tab'

    sj_genes_dict = {}
    for line in open(sj_list, 'r'):
        line = line.strip().split('\t')
        # Fix chromosome
        line[0] = line[0].strip('chr')
        if line[0] == 'X':
            line[0] = '23'
        elif line[0] == 'Y':
            line[0] = '24'
        elif line[0][0] in {'M', 'm'}:
            line[0] = '25'
        elif not re.match(r'^[0-9]+$', line[0]):
            continue
        # Fix sign
        line[3] = '1' if line[3] == '+' else '2'
        idx = tuple(map(int, line[:4]))
        sj_genes_dict[idx] = set(map(int, line[4].split(',')))

    return(sj_genes_dict)

################################################################################
def get_args():
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-i')
    parser.add_argument('-o')
    parser.add_argument('-gd')
    parser.add_argument('-annots')

    return(parser.parse_args())

def main():
    args = get_args()

    annots_dict = parse_annots(args.annots)

    ## Load gene counts
    all_counts_gene_orig, mtx_fns_gene = mmu.load_all_gene_counts(args.i)

    ## Load SJ counts
    sj_annots = load_sj_genes_num(args.gd)
    all_sjs, sj_lists = mmu.get_all_sjs(args.i, sj_annots)
    all_counts_sj, mtx_fns_sj = mmu.load_all_sj_counts(args.i, all_sjs,
        sj_lists)

    ## Make sure both mats are in CSR format
    all_counts_sj = all_counts_sj.tocsr()
    all_counts_gene_orig = all_counts_gene_orig.tocsr()

    ## Need to add on to the mats so there's 1:1 sj:gene
    ##  (some SJs have multiple genes)
    all_counts_sj, all_counts_gene, sj_ord, g_idx = mmu.adjust_mats(
        all_counts_sj, all_counts_gene_orig, all_sjs, sj_annots)

    ## Split up gene and SJ counts
    genes_dict, _ = mmu.split_counts_mat(all_counts_gene, mtx_fns_gene,
        annots_dict)
    sj_dict, _ = mmu.split_counts_mat(all_counts_sj, mtx_fns_sj, annots_dict)
    genes_dict_orig, cell_bcs = mmu.split_counts_mat(all_counts_gene_orig, mtx_fns_gene,
        annots_dict)

    pkl.dump([sj_dict, genes_dict, sj_ord, genes_dict_orig, g_idx, cell_bcs],
        open(args.o, 'wb'))

if __name__ == '__main__':
    main()