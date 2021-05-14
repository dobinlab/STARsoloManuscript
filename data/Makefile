all: whitelists/10Xv2 whitelists/10Xv3

whitelists/10Xv2:
	mkdir -p $(@D) && cd $(@D) && \
	wget http://labshare.cshl.edu/shares/dobin/dobin/STARsolo/Preprint/Data/passlists/10Xv2.gz && \
	gunzip 10Xv2.gz

whitelists/10Xv3:
	mkdir -p $(@D) && cd $(@D) && \
	wget http://labshare.cshl.edu/shares/dobin/dobin/STARsolo/Preprint/Data/passlists/10Xv3.gz && \
	gunzip 10Xv3.gz
