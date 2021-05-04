function genesFNPmean=funPlotGeneDetect(Mas1, Mas, countThreshold, vColor1, vLine1, vH1, vH2, lineW)

%%
m0=full(Mas1)>=countThreshold;
n0=sum(m0,1);
%%
myFigure(2100)
%myFigure(2110)

% myFigure(2200)

for ic=1:length(Mas)

    m1=full(Mas{ic})>=countThreshold;

    n10 = sum( m1 & ~m0, 1);
    n01 = sum(~m1 &  m0, 1);
    n1=sum(m1,1);

    genesFNPmean(ic,1) = mean(n01);
    genesFNPmean(ic,2) = mean(n10);
    genesFNPmean(ic,3) = mean(n01./n1);    
    genesFNPmean(ic,4) = mean(n10./n1);
    
    genesFNPtot(ic,1) = nnz(any(m1,2) & ~any(m0,2)); %nnz(any(m1 & ~m0,2));
     

    %%
    figure(2100)
    %vH=(-25:50:max(n10));

    hh=histogram(n10,vH1, 'FaceColor', vColor1{ic}, 'FaceAlpha', 0.1, 'LineStyle', 'none');
    set(get(get(hh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); %turn off legend for the 2nd line        
    hold on

    [vC]=histcounts(n10,vH1);        
    vHc=vH1(2:end)-(vH1(2)-vH1(1))/2;
    line(vHc,vC,'LineWidth',lineW, 'Color', vColor1{ic}, 'LineStyle', vLine1{ic}) %, 'Marker', 's')        
    box on; grid on
    
    %%
%     figure(2110)
%     %vH=(-25:50:max(n10));
%     vH1p=-0.025:0.05:1.025;
%     
%     hh=histogram(n10./n1, vH1p, 'FaceColor', vColor1{ic}, 'FaceAlpha', 0.1, 'LineStyle', 'none');
%     set(get(get(hh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); %turn off legend for the 2nd line        
%     hold on
% 
%     [vC]=histcounts(n10./n1, vH1p);        
%     vHc=vH1p(2:end)-(vH1p(2)-vH1p(1))/2;
%     line(vHc,vC,'LineWidth',lineW, 'Color', vColor1{ic}, 'LineStyle', vLine1{ic}) %, 'Marker', 's')        
%     box on; grid on

end