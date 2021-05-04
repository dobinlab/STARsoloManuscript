# awk -v cbLen=16 -v geneFile=features.tsv -v whiteList=barcodes.tsv -f truthFromCBUMI.awk

BEGIN {
    ii=0;
    while (getline<geneFile) {
        ii++;
        G[$1]=ii;
        
    };

    ii=0;
    while (getline<whiteList) {
        ii++;
        CB[$1]=ii;
    };

    #print "%\n%";
    #print length(G), length(CB), 0;
}

{
split($0,a,"_");
g=a[2];
if (g=="")
    next; # not genic

getline;
cb=substr($0, 1, cbLen);
umi=substr($0, cbLen+1);

C[G[g]][CB[cb]][umi]=0;

}

END {
for (g in C) {
    for (cb in C[g])
        print g,cb,length(C[g][cb]);
}
}
