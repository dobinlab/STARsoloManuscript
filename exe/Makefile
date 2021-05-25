SHELL=/bin/bash

all: salmon_1.4.0 STAR_2.7.9a generateDecoyTranscriptome.sh mashmap bedtools bwa alevin-fry_0.1.0 kbpy_0.25.0 cellranger-5.0.1

include ../Mf.pars

salmon_%:
	mkdir -p $@_distr && cd $@_distr && \
	wget https://github.com/COMBINE-lab/salmon/releases/download/v$*/salmon-$*_linux_x86_64.tar.gz && \
	tar xvfz salmon-$*_linux_x86_64.tar.gz && \
	cd ../ && ln -s `find ./salmon_$*_distr/ -name salmon` $@

bwa:
	wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 &&\
	tar xvfj bwa-0.7.17.tar.bz2 && 'rm' -f bwa-0.7.17.tar.bz2 && \
	cd bwa-0.7.17 && make && cd .. &&\
	ln -s bwa-0.7.17/bwa bwa

STAR_%:
	curl -o STAR_$* https://raw.githubusercontent.com/alexdobin/STAR/$*/bin/Linux_x86_64_static/STAR && \
	chmod u+x STAR_$*

kbpy_0.25.0:
	conda create -y -n kbpy_0.25.0 python==3.8 pip 
	wget https://github.com/pachterlab/kb_python/archive/refs/tags/v0.25.0.tar.gz &&\
	tar xvfz v0.25.0.tar.gz && rm v0.25.0.tar.gz
	eval "$$(conda shell.bash hook)" && conda activate kbpy_0.25.0 &&\
	cd kb_python-0.25.0 &&\
	python setup.py install

generateDecoyTranscriptome.sh:
	wget https://github.com/COMBINE-lab/SalmonTools/raw/master/scripts/generateDecoyTranscriptome.sh &&\
	chmod u+x generateDecoyTranscriptome.sh
	
mashmap:
	wget https://github.com/marbl/MashMap/releases/download/v2.0/mashmap-Linux64-v2.0.tar.gz &&\
	tar xvfz mashmap-Linux64-v2.0.tar.gz && 'rm' -f mashmap-Linux64-v2.0.tar.gz &&\
	ln -s mashmap-Linux64-v2.0/mashmap
	

bedtools:
	wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary &&\
	chmod u+x bedtools.static.binary &&\
	ln -s bedtools.static.binary bedtools

alevin-fry_0.1.0:
	conda create -y -n alevin-fry_0.1.0 python==3.8 pip
	eval "$$(conda shell.bash hook)" && conda activate alevin-fry_0.1.0 && \
	conda install -y -c bioconda alevin-fry=0.1.0 &&\
	mkdir -p alevin-fry_0.1.0

cellranger-5.0.1:
	wget http://labshare.cshl.edu/shares/dobin/dobin/STARsolo/Preprint/10X/cellranger-5.0.1.tar.gz &&\
	tar xvfz cellranger-5.0.1.tar.gz && 'rm' -f cellranger-5.0.1.tar.gz

gffread:
	wget https://github.com/gpertea/gffread/releases/download/v0.12.6/gffread-0.12.6.Linux_x86_64.tar.gz &&\
	tar xvfz gffread-0.12.6.Linux_x86_64.tar.gz &&\
	'rm' -f gffread-0.12.6.Linux_x86_64.tar.gz &&\
	ln -s gffread-0.12.6.Linux_x86_64/gffread .
