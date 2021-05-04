%function funPlotGeneDetect(Mas1, Mas, countThreshold, vColor1, vLine1)

%%
cellThreshold = 20;
clear n10 n01

n10 = zeros(length(Mas), length(Mas));
n01 = n10;

for ic0=1:length(Mas)
    fprintf(1, '%i ', ic0);
    
    m0=sum(full(Mas{ic0})>=countThreshold,2);

    for ic=ic0+1:length(Mas)

        m1=sum(full(Mas{ic})>=countThreshold,2);

        n10(ic0,ic) = nnz( m1>=cellThreshold & m0==0);
        n10(ic,ic0) = n10(ic0,ic); 
        n01(ic0,ic) = sum( m0>=cellThreshold & m1==0);
        n01(ic,ic0) = n01(ic0,ic); 

    end
end