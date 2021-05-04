function mCorr=funCorrBetweenCellsPlot(mCorr, vLine1, vColor1, lineW)

%%
icorr=2; % Spearman
myFigure(500+icorr)
for ic=1:size(mCorr,3)

    vH=(0:0.01:1);

    hh=histogram(mCorr(:,icorr,ic),vH, 'FaceColor', vColor1{ic}, 'FaceAlpha', 0.1, 'LineStyle', 'none');
    set(get(get(hh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); %turn off legend for the 2nd line        
    hold on

    [vC]=histcounts(mCorr(:,icorr,ic),vH);        
    vHc=vH(2:end)-(vH(2)-vH(1))/2;
    line(vHc,vC,'LineWidth',lineW, 'Color', vColor1{ic}, 'LineStyle', vLine1{ic}) %, 'Marker', 's')        

end
box on; grid on;
