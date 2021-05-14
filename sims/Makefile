include ../Mf.common
.SECONDARY:   #keep all intermediate files

all:
	$(MAKE) sim_mult   genome=human_CR_3.0.0          sample=10X/3/pbmc_5k             cbumiLen=28   cbLen=16   readLen=91   WLfile=10Xv3

extras:
	$(MAKE) sim_mult   genome=h38_Gencode37           sample=10X/3/pbmc_5k             cbumiLen=28   cbLen=16   readLen=91   WLfile=10Xv3
	$(MAKE) sim_mult   genome=mm10_GencodeM21_Diane   sample=10X/2/mouse_ENCSR642XVO   cbumiLen=26   cbLen=16   readLen=151  WLfile=10Xv2

############################################################################################ sims pipeline
.PHONY: sim_mult
sim_mult: \
          ../genomes/$(genome)/genome_transcripts.fa \
          ../genomes/$(genome)/genome_transcripts.1line.fa \
          ../genomes/index/bwa/$(genome)/trgen/BWA_INDEX \
          fastq/$(sample)/cDNA_CBUMI.fq \
          bwa_map/$(genome)/$(sample)/BWA_MAPPED_UNMAPPED \
          bwa_map/$(genome)/$(sample)/MultiGeneNo/SIMS_SAMPLE \
          bwa_map/$(genome)/$(sample)/MultiGeneYes/SIMS_SAMPLE

          #bwa_map/$(genome)/$(sample)/BWA_MAP \
          #bwa_map/$(genome)/$(sample)/MAPPED_UNMAPPED \


%/SIMS_SAMPLE: dirR1R2=$(dd)/samples/$(sample)_sims_$(genome)_$(*F)
%/SIMS_SAMPLE: %/SIMS_FASTQS \
               %/SIMS_R1R2
	$(MAKE) $(dirR1R2)/truth.mtx &&\
	$(MAKE) $(dirR1R2)_OnlyExonicReads


%/truth.mtx:
	cd $* &&\
	awk -O -b -v cbLen=$(cbLen) -v geneFile=$(dd)/genomes/$(genome)/genes.tsv -v whiteList=$(dd)/data/whitelists/$(WLfile) \
		-f $(dd)/sims/truthFromCBUMI.awk _R1_.fq | sort -k2,2n -k1,1n > truth.mtx

%_OnlyExonicReads:
	mkdir -p $@ &&\
	'cp' -f $*/truth.mtx $@ &&\
	grep --no-group-separator ENS -A3 $*/_R2_.fq > $@/_R2_.fq &&\
	grep --no-group-separator ENS -A3 $*/_R1_.fq > $@/_R1_.fq


# generate FASTQs from bwa 
%/MultiGeneNo/SIMS_FASTQS: multOpt=no
%/MultiGeneYes/SIMS_FASTQS: multOpt=yes
%/SIMS_FASTQS:
	mkdir -p $* && cd $* && $(prepComm) &&\
	awk -O -b -v errRate=0.005 -v readLen=$(readLen) -v bcLen=$(cbumiLen) -v multiGene=$(multOpt) -f $(dd)/sims/CBUMIcDNAfqFromMapped_GenicGenomic.awk \
		$(dd)/genomes/$(genome)/genome_transcripts.1line.fa \
		$(dd)/genomes/$(genome)/transcript_to_gene.txt \
		$(dd)/genomes/$(genome)/genes.tsv \
		$(dd)/data/whitelists/$(WLfile) \
		../Mapped.txt >& log.SIMS_FASTQ &&\
	echo "Done" > SIMS_FASTQS

# link FASTQs to R1 R2 in samples/
%/SIMS_R1R2:
	mkdir -p $(dirR1R2) && cd $(dirR1R2) &&\
	ln -fs $(dd)/sims/$*/cDNA.fq  _R2_.fq &&\
	ln -fs $(dd)/sims/$*/CBUMI.fq _R1_.fq

# record CBUMI read from real data into readID of cDNA read for bwa
fastq/%/cDNA_CBUMI.fq: sample=$*
fastq/%/cDNA_CBUMI.fq:
	mkdir -p $(@D) && cd $(@D) && \
	paste  <(zcat $(dd)/samples/$(sample)/*_R1_*) <(zcat $(dd)/samples/$(sample)/*_R2_*) | \
	awk -O -b '{getline; cbumi=substr($$1,1,$(cbumiLen)); cdna=$$2; getline;getline; ii++; name= "@" ii "_" cbumi;  print name "\n" cdna "\n+\n" $$2 > "cDNA_CBUMI.fq"}' && \
	echo "Done" > CO

# combine genome.fa and transcripts.fa, and make them 1line fasta
%/genome_transcripts.1line.fa:
	cd $* &&\
        cat  genome.fa transcripts.fa | awk -O -b '{if ($$1~/^>/) $$0="\n" substr($$1,2) " "; printf $$0}' > genome_transcripts.1line.fa
%/genome_transcripts.fa:
	cd $* &&\
        cat  genome.fa transcripts.fa > genome_transcripts.fa


# bwa index
%/BWA_INDEX: $(dd)/genomes/$(genome)/genome_transcripts.fa
	mkdir -p $* && cd $* &&\
	$(dd)/exe/bwa index -p bwaIndex $< >& log.index &&\
	echo "Done" > BWA_INDEX


# bwa map to trgen
bwa_map/%/BWA_MAP:
	mkdir -p $(@D) && cd $(@D) &&\
	$(dd)/exe/bwa mem -a -M -t 20    $(dd)/genomes/index/bwa/$(genome)/trgen/bwaIndex    $(dd)/sims/fastq/$(sample)/cDNA_CBUMI.fq    1>aligned.sam 2>log.map    &&\
	echo "Done" > BWA_MAP

# Extract mapped reads from SAM. Also unmapped, but they are not used.
%/MAPPED_UNMAPPED:
	cd $* && $(prepComm) &&\
	cat aligned.sam | awk -O -b -v bcLen=$(cbumiLen) -v cbLen=$(cbLen) -f $(dd)/sims/samToMappedAndUnmapped_mult.awk &&\
	echo Done > MAPPED_UNMAPPED

# bwa map and convert: save disk space. Replaces BWA_MAP and MAPPED_UNMAPPED
bwa_map/%/BWA_MAPPED_UNMAPPED:
	mkdir -p $(@D) && cd $(@D) &&\
	$(dd)/exe/bwa mem -a -M -t 20    $(dd)/genomes/index/bwa/$(genome)/trgen/bwaIndex    $(dd)/sims/fastq/$(sample)/cDNA_CBUMI.fq  2>log.map \
		| awk -O -b -v bcLen=$(cbumiLen) -v cbLen=$(cbLen) -f $(dd)/sims/samToMappedAndUnmapped_mult.awk   &&\
	echo "Done" > BWA_MAPPED_UNMAPPED



