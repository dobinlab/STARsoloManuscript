all: pbmc5k/wallTime.txt pbmc5k/residentRAM.txt pbmc5k_pre/wallTime.txt pbmc5k_pre/residentRAM.txt

%/wallTime.txt: grepText='Elapsed (wall clock) time'
%/residentRAM.txt: grepText='Maximum resident set size (kbytes)'

pbmc5k/%:
	mkdir -p $(@D)
	# STAR
	grep $(grepText) ../count/STAR_2.7.9x/human_CR_3.0.0/*/10X_CR4_noSAM/10X/3/pbmc_5k/*/b*/log | awk '{split($$1,a,"/"); print a[3],a[5],a[6],a[10],a[11],"-", $$NF}' > $@
        # kb
	grep $(grepText) ../count/kbpy_0.25.0/human_CR_3.0.0/standard_1/default/10X/3/pbmc_5k/*/b*/log | awk '{split($$1,a,"/"); print a[3],a[5],a[6],a[10],a[11],"-", $$NF}' >> $@
	grep $(grepText) ../count/CellRanger_5.0.1/human_CR_3.0.0/standard/default/10X/3/pbmc_5k/*/b*/log | awk '{split($$1,a,"/"); print a[3],a[5],a[6],a[10],a[11],"-", $$NF}' >> $@
	# salmon
	grep $(grepText) ../count/*/human_CR_3.0.0/*/{rad,sketch_rad}/10X/3/pbmc_5k/*/b*/log.salmon | awk '{split($$1,a,"/"); print a[3],a[5],a[6],a[10],a[11],"-", $$NF}' >> $@
	# alevin counts
	grep $(grepText) ../count/*/human_CR_3.0.0/*/{rad,sketch_rad}/10X/3/pbmc_5k/*/b*/gpl_knee/log.{gpl,collate,quant_cr-like} | awk '{split($$1,a,"/"); print a[3],a[5],a[6],a[10],a[11],a[13], $$NF}' >> $@

pbmc5k_pre/%:
	mkdir -p $(@D)
	# STAR
	grep $(grepText) ../count/STAR_2.7.9x/human_CR_3.0.0/*/10X_CR4_GeneFull_noSAM/10X/3/pbmc_5k/*/b*/log | awk '{split($$1,a,"/"); print a[3],a[5],a[6],a[10],a[11],"-", $$NF}' > $@
        # kb
	grep $(grepText) ../count/kbpy_0.25.0/human_CR_3.0.0/nucleus_1/default/10X/3/pbmc_5k/*/b*/log | awk '{split($$1,a,"/"); print a[3],a[5],a[6],a[10],a[11],"-", $$NF}' >> $@
