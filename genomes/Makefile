include ../Mf.common
.SECONDARY:   # keep all intermediate files

all: human_CR_3.0.0/GENOME \
     human_CR_3.0.0/INDEXES

h38_Gencode37: \
     h38_Gencode37/GENOME \
     h38_Gencode37/INDEXES \

mouse_CR3: \
     mouse_CR_3.0.0/GENOME \
     mouse_CR_3.0.0/INDEXES

forDiane: \
     mm10_GencodeM21_Diane/GENOME \
     mm10_GencodeM21_Diane/INDEXES

#################################################################################################################################### genome specifics
human_CR_3.0.0/%: CRgenomeDir=$(dd)/genomes/human_CR_3.0.0/refdata-cellranger-GRCh38-3.0.0

human_CR_3.0.0/refdata-cellranger-GRCh38-3.0.0: 
	mkdir -p $(@D) && cd $(@D) && \
	wget -r http://labshare.cshl.edu/shares/dobin/dobin/STARsolo/Preprint/Genomes/refdata-cellranger-GRCh38-3.0.0 && \
	mv labshare.cshl.edu/shares/dobin/dobin/STARsolo/Preprint/Genomes/refdata-cellranger-GRCh38-3.0.0/ . && \
   	rm -rf labshare.cshl.edu

human_CR_3.0.0/genome.fa: human_CR_3.0.0/refdata-cellranger-GRCh38-3.0.0
	mkdir -p $(@D) && cd $(@D) && 'rm' -f genome.fa &&\
	ln -s refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa genome.fa

human_CR_3.0.0/annotations.gtf: human_CR_3.0.0/refdata-cellranger-GRCh38-3.0.0
	mkdir -p $(@D) && cd $(@D) && 'rm' -f annotations.gtf &&\
	ln -s refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf annotations.gtf

############################################
mouse_CR_3.0.0/%: CRgenomeDir=/scratch/dobin/genomes/Mouse/refdata-cellranger-mm10-3.0.0/

mouse_CR_3.0.0/genome.fa:
	mkdir -p $(@D) && cd $(@D) && \
	ln -s /scratch/dobin/genomes/Mouse/refdata-cellranger-mm10-3.0.0/fasta/genome.fa

mouse_CR_3.0.0/annotations.gtf:
	mkdir -p $(@D) && cd $(@D) && \
	ln -s /scratch/dobin/genomes/Mouse/refdata-cellranger-mm10-3.0.0/genes/genes.gtf annotations.gtf

############################################
h38_Gencode37/genome.fa:
	mkdir -p h38_Gencode37 && cd h38_Gencode37 &&\
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz &&\
	gunzip GRCh38.primary_assembly.genome.fa.gz &&\
	ln -s GRCh38.primary_assembly.genome.fa genome.fa

h38_Gencode37/annotations.gtf:
	mkdir -p h38_Gencode37 && cd h38_Gencode37 &&\
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.primary_assembly.annotation.gtf.gz &&\
	gunzip gencode.v37.primary_assembly.annotation.gtf.gz &&\
	ln -s gencode.v37.primary_assembly.annotation.gtf annotations.gtf
############################################
mm10_GencodeM21_Diane/annotations.gtf: mm10_GencodeM21_Diane/mm10-M21_minimal_det-male.gtf
	# add spike-in transcripts to satisfy kb: copy exon lines replacing "exon" with "transcript"
	mkdir -p $(@D) && cd $(@D) &&\
	awk 'BEGIN {FS=OFS="\t"} {print; if ($$2=="spikein") {$$3="transcript"; print} }' mm10-M21_minimal_det-male.gtf > annotations.gtf

mm10_GencodeM21_Diane/annotations.gtf.gz:
	mkdir -p $(@D) && cd $(@D) &&\
	wget http://woldlab.caltech.edu/~diane/encode-202006-jamboree-detrout-rna-sc-pipeline/genome/mm10-M21_minimal-male/mm10-M21_minimal_det-male.gtf

mm10_GencodeM21_Diane/genome.fa: mm10_GencodeM21_Diane/mm10_no_alt_analysis_set_ENCODE.fasta
	mkdir -p $(@D) && cd $(@D) &&\
	wget https://www.encodeproject.org/files/ENCFF335FFV/@@download/ENCFF335FFV.fasta.gz &&\
	gunzip ENCFF335FFV.fasta.gz &&\
	wget https://www.encodeproject.org/files/ENCFF001RTP/@@download/ENCFF001RTP.fasta.gz &&\
	gunzip ENCFF001RTP.fasta.gz &&\
	cat mm10_no_alt_analysis_set_ENCODE.fasta ENCFF001RTP.fasta ENCFF335FFV.fasta > genome.fa

mm10_GencodeM21_Diane/mm10_no_alt_analysis_set_ENCODE.fasta:
	mkdir -p $(@D) && cd $(@D) &&\
	wget https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz &&\
	gunzip mm10_no_alt_analysis_set_ENCODE.fasta.gz

############################################################################################################################################
%/GENOME: %/genome.fa %/annotations.gtf %/transcripts.fa %/transcript_to_gene.txt %/genes.tsv
	echo Done > $(@D)/GENOME

%/transcripts.fa: %/genome.fa %/annotations.gtf
	cd $* && \
	samtools faidx genome.fa && \
	$(dd)/exe/gffread -w transcripts.fa -g genome.fa annotations.gtf

%/transcript_to_gene.txt: %/annotations.gtf
	cd $* && \
	awk '$$3=="transcript" {gene=$$0; gsub(/.*gene_id "/,"",gene); gsub(/".*/,"",gene); tr=$$0; gsub(/.*transcript_id "/,"",tr); gsub(/".*/,"",tr); print tr "\t" gene }' \
	    annotations.gtf > transcript_to_gene.txt && \
	awk '{print $$1 "\t" $$2}' transcript_to_gene.txt > transcript_to_gene.2col.txt

%/genes.tsv: %/transcript_to_gene.txt
	cd $* && \
	awk '{if (!($$2 in G)) print $$2; G[$$2]=0}'  transcript_to_gene.txt > genes.tsv

.PHONY: %/INDEXES
%/INDEXES: genome = $*
%/INDEXES:
	[[ ! -z "$(CRgenomeDir)" ]] && mkdir -p index/CellRanger_5.0.1/$(genome) && cd index/CellRanger_5.0.1/$(genome) && ln --force -s $(CRgenomeDir) standard
	$(MAKE) -f Mf.index index/kbpy_0.25.0/$(genome)/standard_1/CO
	$(MAKE) -f Mf.index index/salmon_1.4.0/$(genome)/standard/CO
	$(MAKE) -f Mf.index index/salmon_1.4.0/$(genome)/decoyFull/CO
	$(MAKE) -f Mf.index index/STAR_2.7.9a/$(genome)/fullSA/CO
	$(MAKE) -f Mf.index index/STAR_2.7.9a/$(genome)/sparseSA3/CO
	$(MAKE) -f Mf.index index/salmon_1.4.0/$(genome)/decoyPartial/CO
	$(MAKE) -f Mf.index index/kbpy_0.25.0/$(genome)/nucleus_1/CO

