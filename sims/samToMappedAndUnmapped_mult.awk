# -v bcLen=28 -v cbLen=16

BEGIN {
    OFS="\t";
    for (ii=1;ii<=bcLen;ii++)  bcQual=bcQual "N";

    scoreField=14;
}

{
    if (substr($1,1,1)=="@") {
        print > "Unmapped.sam";
        next;
    };

    if ($3=="*") {
        print > "Unmapped.sam";
        print "@" $1 "\n" $10 "\n+\n" $11 > "Unmapped_cDNA.fq";
        print "@" $1 "\n" cbumi "\n+\n" bcQual > "Unmapped_CBUMI.fq";
        next;
    };

    score1=substr($scoreField,6);

    if (r=="") {# very first read
        r=$1;
        scoreTop=score1;
    };

    if (r!=$1) {
        cbumi=substr(r,length(r)-bcLen+1);
        #if (index(cbumi,"N")>0) next;
        print substr(cbumi,1,cbLen), substr(cbumi,cbLen+1), "G", substr(tt,2), r, score1 > "Mapped.txt";

        tt="";
        r=$1;
        scoreTop=score1; # bwa output top score read first
    };
    
    if (scoreTop==score1) { # keep all reads with score=scorTop
        tt=tt ";" $3 "," $4 ",0";
    };
}
END {
    cbumi=substr(r,length(r)-bcLen+1);
    print substr(cbumi,1,16), substr(cbumi,17), "G", substr(tt,2), r, score1 > "Mapped.txt";
}
