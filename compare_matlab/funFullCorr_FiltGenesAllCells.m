function corrM=funFullCorr_FiltGenesAllCells(Ma, countThr, nCorr, label1)

nCases=length(Ma);
if nCorr==0
    nCorr=nCases;
end

corrM=ones(nCases, nCases, 2);

for ic=1 : nCorr
    m=full(Ma{ic}(:));
    if isempty(m)
        continue;
    end
    %%
    vv = m >= countThr ;
    for ic1=(ic+1):nCases
        fprintf(1, '%i ', ic1);

        m1=full(Ma{ic1}(:));
        if isempty(m1)
            continue;
        end
        
        vv1 = vv | m1 >= countThr;

        
        %corrM(ic,ic1,1) = corr(log(m+1), log(m1+1));
        corrM(ic,ic1,2) = corr(m(vv1), m1(vv1), 'type','Spearman');
        
        corrM(ic1,ic,:) =corrM(ic,ic1,:);
        
        disp([label1{ic} '   ' label1{ic1} '   ' num2str(corrM(ic,ic1,2))])
        
    end
end
