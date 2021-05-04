function vCorr=funCorrBetweenCells_FiltGenesAllCells(x,y, countMin)
% computer correlation between columns in 2 matrices
nCB=size(x,2);
vCorr=zeros(nCB,2);

% keep genes that are expressed in at least one cell
vv=any(x>=countMin,2) | any(y>=countMin,2);

%vv=sum(x>=countMin,2)>=10 | any(y>=countMin,2)>=10;


for ic=1:nCB
    %if mod(ic,1000)==0; fprintf(1,"%i ",ic); end
    
    x1=x(vv,ic);
    y1=y(vv,ic);
    
    nx=sum(x1);
    ny=sum(y1);

    cc=corr([log2(x1/nx*1e4+1) log2(y1/ny*1e4+1)], 'type', 'Pearson');
    vCorr(ic,1)=cc(1,2);
    cc=corr([x1 y1], 'type', 'Spearman');
    vCorr(ic,2)=cc(1,2);    
end
