function vCorr=funCorrMat(x,y, countMin)
% computer correlation between columns in 2 matrices
nCB=size(x,2);
vCorr=zeros(nCB,2);
for ic=1:nCB
    x1=x(:,ic);
    y1=y(:,ic);
    if mod(ic,1000)==0; fprintf(1,"%i ",ic); end
    vv=x1>=countMin | y1>=countMin;
    x1=x1(vv);
    y1=y1(vv);

    cc=corr([log2(x1+1) log2(y1+1)], 'type', 'Pearson');
    vCorr(ic,1)=cc(1,2);
    cc=corr([x1 y1], 'type', 'Spearman');
    vCorr(ic,2)=cc(1,2);    
end
