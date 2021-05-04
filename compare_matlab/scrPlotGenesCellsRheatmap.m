label1=casesSelectNamesShort;
myFigure(101)
h=heatmap(label1,label1, round(corrGenesCells(vSelect,vSelect,2),3), 'FontSize',24, 'ColorbarVisible','off');
px=0.3; py=0.24; 
set(gca, 'Position', [px py 0.99-px 0.95-py]);

s=struct(h);
s.XAxis.TickLabelRotation = 45;

%title('Spearman R for gene/cell counts')
%s.TitleHandle.FontWeight = 'normal';

saveas(gcf, [savePrefix 'SpearmanHeatmap' '.png']);
