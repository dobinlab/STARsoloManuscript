SHELL=/bin/bash

all: salmon_1.4.0 STAR_2.7.9a kbpy_0.25.0 generateDecoyTranscriptome.sh mashmap bedtools

include ../Mf.pars

salmon_%:
	mkdir -p $@_distr && cd $@_distr && \
	wget https://github.com/COMBINE-lab/salmon/releases/download/v$*/salmon-$*_linux_x86_64.tar.gz && \
	tar xvfz salmon-$*_linux_x86_64.tar.gz && \
	cd ../ && ln -s `find ./salmon_$*_distr/ -name salmon` $@

bwa:
	@echo need to download bwa

STAR_%:
	curl -o STAR_$* https://raw.githubusercontent.com/alexdobin/STAR/2.7.5a/bin/Linux_x86_64_static/STAR && \
	chmod u+x STAR_$*

kbpy_0.25.0:
	# worked from command line!
	conda create -n kbpy_0.25.0 python==3.8 pip
	conda activate kbpy_0.25.0
	python setup.py install

generateDecoyTranscriptome.sh:
	wget https://github.com/COMBINE-lab/SalmonTools/raw/master/scripts/generateDecoyTranscriptome.sh &&\
	chmod u+x generateDecoyTranscriptome.sh
	
mashmap:
	wget https://github.com/marbl/MashMap/releases/download/v2.0/mashmap-Linux64-v2.0.tar.gz &&\
	tar xvfz mashmap-Linux64-v2.0.tar.gz &&\
	ln -s mashmap-Linux64-v2.0/mashmap
	

bedtools:
	wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary &&\
	chmod u+x bedtools.static.binary &&\
	ln -s bedtools.static.binary bedtools
