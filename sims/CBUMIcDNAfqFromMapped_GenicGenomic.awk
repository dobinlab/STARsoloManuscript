# usage awk -v errRate=0.005 -v readLen=91 -v bcLen=28 [-v multiGene=no] -f .awk Transcripts.fa TranscriptGene.txt GeneID.txt CB.txt genic.txt

BEGIN {
    #readLen=98;
    for (ii=1;ii<=bcLen;ii++) qual1 = qual1 "F";
    for (ii=1;ii<=readLen;ii++) qual2 = qual2 "F";
    for (ii=1;ii<=readLen;ii++) polyA = polyA "A";

    errRate+=0;
    errRate*=4/3;

    split("ACGT",nucl,"");

    stderr="/dev/stderr";
    OFS="\t";

    print errRate, readLen, bcLen, multiGene > stderr;
}
BEGINFILE {
    if (ARGIND==5) {
        print "trSeq: " length(trSeq) > stderr;
        print "G: " length(G) > stderr;
        print "G: " length(G) > stderr;
        print "CB: " length(CB) > stderr;
   }
}
{
   if (ARGIND==1) {
       trSeq[$1] = $2 polyA;

   } else if (ARGIND==2) {
       kk++;
       TG[$1]=$2;

   } else if (ARGIND==3) {
       mm++;
       G[$1]=mm;

   } else if (ARGIND==4) {
       jj++;
       CB[$1]=jj;

   } else {

       cb=$1;umi=$2;
       if ( !(cb in CB) || index(umi,"N")>0)
           next; # cb should be present in WL, do not allow N in UMI

       split($4,aligns,";"); # number of mapped loci, tr and genomic

       split("",trs); # this will contain transcript tuples vs index
       split("",trst);# same vs trID
       split("",genes); # genes
       for (ii in aligns) {
           if (aligns[ii]~/^ENS/) {
               split(aligns[ii],b,","); # b[1]=trID
               trs[length(trs)+1]=aligns[ii];
               trst[b[1]]=aligns[ii];
               genes[TG[b[1]]]=0;
           };
       };

       print length(aligns), length(trs), length(genes) > stderr;

       if (length(trs)==0) {# no transcripts
           if (cbUmiTr[cb umi]=="") {# new cb/umi, mark it with non-exonic
               cbUmiTr[cb umi] = "non-exonic";
           } else if (cbUmiTr[cb umi]!="non-exonic") {# this cb/umi was already marked as transcript, skip this read
               print cbUmiTr[cb umi], $0 > "Conflicts.txt";
               next;
           };
           n1 = int(length(aligns)*rand()+1);
           split(aligns[n1],align1s,",");

       } else {# transcripts are present, the align will be chosen from them
           if (multiGene=="no" && length(genes)>1 ) {
               print "nGenes=" length(genes), $0 > "MultiGenicReads.txt";
               next; # multigene read - do not record
           }; 


           if (cbUmiTr[cb umi] =="") {# first time - will assign one transcript to this UMI randomly
               # select randomly one transcript/gene from tr
               n1 = int(length(trs)*rand()+1);
               align1=trs[n1];
               split(align1,align1s,",");
               cbUmiTr[cb umi] = align1s[1];               
           } else if (cbUmiTr[cb umi] =="non-exonic") {# this cb/umi was marked as transcript before, skip this read
               next;
           } else {# find transcript that was already assigned
               align1=trst[cbUmiTr[cb umi]];

               if (align1=="") {
                   print "nTr=" length(trs), trs[1], cbUmiTr[cb umi], $0 > "Conflicts.txt";
                   next;
               };

               split(align1,align1s,",");  
           };
       };

       id = align1s[1]; # id of the align from which read will be simulated: could be tr or chr
       seq1 = substr(trSeq[id],align1s[2]+0,readLen);

       if (errRate>0) {
           for (ii=1;ii<=length(seq1);ii++) {
               if (rand()<errRate) {
                   b1=nucl[int(rand()*4)+1];
                   if (b1!=substr(seq1,ii,1))
                       nErr++;
                   seq1=substr(seq1,1,ii-1) b1 substr(seq1,ii+1);
               };
           };
       };

       nBases+=length(seq1);

       print "@" readNamePrefix FNR "_" TG[id] "_" id "_" b[2] "\n" seq1 "\n+\n" qual2 > "cDNA.fq";
       print "@" readNamePrefix FNR "_" TG[id] "_" id "\n" cb umi "\n+\n" qual1 > "CBUMI.fq";

       if (length(trs)>0) 
           nUMI[G[TG[id]]][CB[cb]][umi]=0;
   };
};

END {
    print nBases, nErr, nErr/nBases > stderr;
    for (gg in nUMI) {
        for (cb in nUMI[gg]) {
            print gg,cb, length(nUMI[gg][cb]) > "trueUMIcounts.mtx";
        };
    };
};
