# collect Elapsed time
# usage: 

BEGIN {
    OFS=","
}

{
sample=$1; # sample
tin=$NF; # time

gsub("/b[0-9]/","",sample);

if (sample!=sampleOld){
    if (sampleOld!="")
        print sampleOld, N, St/N, (N<2 ? -1 : ((St2-St^2/N)/(N-1))^(1/2)) tAll;
    N=St=St2=0;
    tAll="";
    sampleOld=sample;
};

n=split(tin,tt,":");
t=tt[n];
t+=tt[n-1]*60;
t+=tt[n-2]*3660;

St+=t;
St2+=t^2;
N++;
tAll=tAll OFS t;

#print sample,St,St2,N,tAll;

}
END {
    print sampleOld, N, St/N, (N==0 ? -1 : ((St2-St^2/N)/(N-1))^(1/2)) tAll;
}
