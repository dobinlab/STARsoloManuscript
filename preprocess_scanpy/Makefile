alevin/%:
	cd $* && \
	cp /scratch/dobin/STAR/STARsoloPreprint/count/CellRanger_5.0.1/human_CR_3.0.0/standard/default/10X/3/pbmc_5k/20/b02/Run1/outs/filtered_feature_bc_matrix/features.tsv.gz .
	
	nGenes=`zcat features.tsv.gz | wc -l`
	awk -v nGenes=$nGenes '{if (NR==3) $2=nGenes; if (NR>2 && $2<=nGenes) {f1=$1; $1=$2; $2=f1}; print}' quants_mat.mtx | gzip -c > matrix.mtx.gz
	gzip -c quants_mat_rows.txt > barcodes.tsv.gz


kb:
	cp /scratch/dobin/STAR/STARsoloPreprint/count/CellRanger_5.0.1/human_CR_3.0.0/standard/default/10X/3/pbmc_5k/20/b02/Run1/outs/filtered_feature_bc_matrix/features.tsv.gz .
	awk '{if (NR>3) {f1=$1; $1=$2; $2=f1}; print}' cells_x_genes.mtx | gzip -c > matrix.mtx.gz
	gzip -c cells_x_genes.barcodes.txt > barcodes.tsv.gz 

