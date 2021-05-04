function [RD,MARD]=funCalcRelDiff(Ma,countThr)

%
for ic=1   %:length(M) % ic=1 ref-tool to compare with
    m=full(Ma{ic}(:));
    if isempty(m)
        continue;
    end
    vv = m >= countThr;
    
    for ic1=ic+1:length(Ma)

        fprintf(1, '%i ', ic1);

        m1=full(Ma{ic1}(:));
        if isempty(m1)
            continue;
        end
        
        vv1 = vv | m1 >= countThr;
        m2=m(vv1);
        m1=m1(vv1);        
        
        ma1=(m1-m2)./max(m2,m1);
        
        RD{ic1,ic}=ma1;
        MARD(ic1,ic)=mean(abs(ma1));
    end
end