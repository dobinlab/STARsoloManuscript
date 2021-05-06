include Mf.common

.SECONDARY:

.PHONY: all

all: prep
#all: prep real sims bench bench_memUsage



real:
	$(MAKE) All_Gene_Real       genome=human_CR_3.0.0   sample1=10X/3/pbmc_5k   threadRun=20/run1
	$(MAKE) All_GeneFull_Real   genome=human_CR_3.0.0   sample1=10X/3/pbmc_5k   threadRun=20/run1

sims:
	$(MAKE) All_Gene_Real       genome=human_CR_3.0.0   sample1=10X/3/pbmc_5k_sims_MultiGeneNo    threadRun=20/run1
	$(MAKE) All_Gene_Real       genome=human_CR_3.0.0   sample1=10X/3/pbmc_5k_sims_MultiGeneYes   threadRun=20/run1
	$(MAKE) All_Gene_Real       genome=human_CR_3.0.0   sample1=10X/3/pbmc_5k_sims_MultiGeneNo_OnlyExonicReads   threadRun=20/run1
	$(MAKE) All_Gene_Real       genome=human_CR_3.0.0   sample1=10X/3/pbmc_5k_sims_MultiGeneYes_OnlyExonicReads  threadRun=20/run1

bench:
	$(MAKE) bench/All_Gene_Real benchThreadNumbers="20 8 12 4 16" genome=human_CR_3.0.0 sample=10X/3/pbmc_5k
	$(MAKE) bench/All_GeneFull_Real benchThreadNumbers="20 8 12 4 16" genome=human_CR_3.0.0 sample=10X/3/pbmc_5k 


bench_memUsage:
	$(MAKE) bench/Gene/CellRanger benchThreadNumbers="20 16 12 8 4" benchRunIndexes="memUse1" genome=human_CR_3.0.0 sample=10X/3/pbmc_5k PREPoption="memUsage"


####################################################################################################################### extras
real_mult:
	$(MAKE) -f Mf.count count/STAR_2.7.9x/human_CR_3.0.0/fullSA/10X_CR4_noSAM_mult/10X/3/pbmc_5k/20/run1/CO
	$(MAKE) -f Mf.count count/STAR_2.7.9x/human_CR_3.0.0/fullSA/10X_CR4_noSAM_mult_UniqueGenomic/10X/3/pbmc_5k/20/run1/CO

release_tests:
	$(MAKE) -f Mf.count count/STAR_2.7.9a/human_CR_3.0.0/fullSA/10X_CR4_noSAM/10X/3/pbmc_5k/20/relTest/CO
	$(MAKE) -f Mf.count count/STAR_2.7.9a/human_CR_3.0.0/fullSA/10X_CR4_noSAM_mult/10X/3/pbmc_5k/20/relTest/CO
	$(MAKE) -f Mf.count count/STAR_2.7.9a/human_CR_3.0.0/fullSA/10X_noSAM_sims_mult_ENCODE/10X/3/pbmc_5k_sims_MultiGeneNo/20/relTest/CO
	$(MAKE) -f Mf.count count/STAR_2.7.9a/human_CR_3.0.0/fullSA/10X_noSAM_sims_mult_ENCODE/10X/3/pbmc_5k_sims_MultiGeneYes/20/relTest/CO

real_TM_select5: samples=TM_10X_P7_0 TM_10X_P4_2 TM_10X_P7_9 TM_10X_P8_14 TM_10X_P7_6
real_TM_select5:
	for s1 in $(samples); do \
		$(MAKE) Gene_Real_Ssp-Afd-Ask-kb genome=mouse_CR_3.0.0 sample=10X/2/$$s1 threadRun=20/run1; \
	done

sims_all_mouse:
	$(MAKE) All_Gene_Sims   genome=mm10_GencodeM21_Diane   sample=10X/2/mouse_ENCSR642XVO_sims_mm10_GencodeM21_Diane_MultiGeneNo                    threadRun=20/run1
	$(MAKE) All_Gene_Sims   genome=mm10_GencodeM21_Diane   sample=10X/2/mouse_ENCSR642XVO_sims_mm10_GencodeM21_Diane_MultiGeneYes                   threadRun=20/run1
	$(MAKE) All_Gene_Sims   genome=mm10_GencodeM21_Diane   sample=10X/2/mouse_ENCSR642XVO_sims_mm10_GencodeM21_Diane_MultiGeneNo_OnlyExonicReads    threadRun=20/run1
	$(MAKE) All_Gene_Sims   genome=mm10_GencodeM21_Diane   sample=10X/2/mouse_ENCSR642XVO_sims_mm10_GencodeM21_Diane_MultiGeneYes_OnlyExonicReads   threadRun=20/run1


############################################################ Versions
alevin-fry := salmon-alevin-fry_1.4.0_0.1.0
kbpy       := kbpy_0.25.0
STAR       := STAR_2.7.9a
CellRanger := CellRanger_5.0.1

############################################################ Bench
#benchThreadNumbers := 20 12 8 4 16
benchRunIndexes    ?= b01 b02 b03 b04 b05


bench/%: benchOpts= PREPoption=sleep
bench/%:
	for rr in $(benchRunIndexes); do \
	for tt in $(benchThreadNumbers); do \
		$(MAKE) $* threadRun=$$tt/$$rr genome=$(genome) sample=$(sample) $(benchOpts) ; \
	done; \
	done
		#$(MAKE) $* threadRun=$$tt/$$rr genome=$(genome) sample=$(sample); \

############################################################ tools and options
# Gene: standard gene expression

All_GeneFull_Real:   \
                     GeneFull/kbpy \
                     GeneFull/STAR_CR4 \
                     # GeneFull/CellRanger \
	
All_Gene_Sims: Gene/STAR_sims \
               Gene/kbpy \
               Gene/alevin-fry

All_Gene_Real:  \
                Gene/kbpy \
                Gene/alevin-fry \
                Gene/STAR_CR4 \
                Gene/CellRanger

Gene_Real_Ssp-Afd-Ask-kb:
	make -f Mf.count count/$(STAR)/$(genome)/sparseSA3/10X_CR4_noSAM/$(sample)/$(threadRun)/CO
	make -f Mf.count count/$(kbpy)/$(genome)/standard_1/default/$(sample)/$(threadRun)/CO
	make -f Mf.count count/$(alevin-fry)/$(genome)/standard/sketch_rad/$(sample)/$(threadRun)/CO afGplOpt=knee
	make -f Mf.count count/$(alevin-fry)/$(genome)/decoyFull/rad/$(sample)/$(threadRun)/CO       afGplOpt=knee

Gene/CellRanger: 
	make -f Mf.count count/$(CellRanger)/$(genome)/standard/default/$(sample)/$(threadRun)/CO

GeneFull/CellRanger: 
	make -f Mf.count count/$(CellRanger)/$(genome)/standard/intronic/$(sample)/$(threadRun)/CO

Gene/STAR_sims:
	$(MAKE) -f Mf.count count/$(STAR)/$(genome)/fullSA/10X_noSAM_sims_mult_ENCODE/$(sample)/$(threadRun)/CO

Gene/STAR_CR4:
	make -f Mf.count count/$(STAR)/$(genome)/fullSA/10X_CR4_noSAM/$(sample)/$(threadRun)/CO
	make -f Mf.count count/$(STAR)/$(genome)/sparseSA3/10X_CR4_noSAM/$(sample)/$(threadRun)/CO

GeneFull/STAR_CR4:
	$(MAKE) -f Mf.count count/$(STAR)/$(genome)/fullSA/10X_CR4_GeneFull_noSAM/$(sample)/$(threadRun)/CO
	$(MAKE) -f Mf.count count/$(STAR)/$(genome)/sparseSA3/10X_CR4_GeneFull_noSAM/$(sample)/$(threadRun)/CO

Gene/kbpy:
	make -f Mf.count count/$(kbpy)/$(genome)/standard_1/default/$(sample)/$(threadRun)/CO
	make -f Mf.count count/$(kbpy)/$(genome)/standard_1/mult/$(sample)/$(threadRun)/CO

GeneFull/kbpy:
	make -f Mf.count count/$(kbpy)/$(genome)/nucleus_1/default/$(sample)/$(threadRun)/CO
	#make -f Mf.count count/$(kbpy)/$(genome)/nucleus_8/default/$(sample)/$(threadRun)/CO
	#make -f Mf.count count/$(kbpy)/$(genome)/nucleus_2/default/$(sample)/$(threadRun)/CO
	#make -f Mf.count count/$(kbpy)/$(genome)/nucleus_4/default/$(sample)/$(threadRun)/CO

Gene/alevin-fry:
	make -f Mf.count count/$(alevin-fry)/$(genome)/standard/sketch_rad/$(sample)/$(threadRun)/CO afGplOpt=knee
	make -f Mf.count count/$(alevin-fry)/$(genome)/standard/rad/$(sample)/$(threadRun)/CO        afGplOpt=knee
	make -f Mf.count count/$(alevin-fry)/$(genome)/decoyFull/rad/$(sample)/$(threadRun)/CO       afGplOpt=knee
	make -f Mf.count count/$(alevin-fry)/$(genome)/decoyPartial/rad/$(sample)/$(threadRun)/CO    afGplOpt=knee


################################## Prep
.PHONY: prep exe genomes samples genomes/index

prep: samples exe

exe:
	$(MAKE) -C exe

genomes:
	$(MAKE) -C genomes

samples:
	$(MAKE) -C samples

