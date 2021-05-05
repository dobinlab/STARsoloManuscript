function funPlotHists(vH, M, vColor1, vLine1, iFig)

%%
myFigure(iFig)
myFigure(2200)

for ic=1:size(M,2)

    figure(iFig)

    hh=histogram(M(:,ic),vH, 'FaceColor', vColor1{ic}, 'FaceAlpha', 0.1, 'LineStyle', 'none');
    set(get(get(hh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); %turn off legend for the 2nd line        
    hold on

    [vC]=histcounts(M(:,ic),vH);        
    vHc=vH(2:end)-(vH(2)-vH(1))/2;
    line(vHc,vC,'LineWidth',4, 'Color', vColor1{ic}, 'LineStyle', vLine1{ic}) %, 'Marker', 's')        
    box on; grid on

end