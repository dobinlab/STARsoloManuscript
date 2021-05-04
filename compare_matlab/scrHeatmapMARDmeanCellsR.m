%% MARD / mean-Spearman-cells table
corrCellsMean = squeeze(mean(corrCells(:,2,vSelect)));

MARDselect = MARD(vSelect);
label1=casesSelectNamesShort;
myFigure(601)
axes
set(gca, 'Position', [0.45 0.1 0.25 0.8]);
h=heatmap({'\bf MARD'}, label1, round([MARDselect],3), 'FontSize',32, 'ColorbarVisible','off');
axes
s=struct(h);
s.Axes.XAxisLocation = 'top';
h.Colormap = h.Colormap(end:-1:1,:);
s.XAxis.TickLabelInterpreter = 'latex';
% s.YAxis.FontSize=28;

set(gca, 'Position', [0.71 0.1 0.25 0.8]);
h=heatmap({'\bf Mean Cells R'}, label1, round([corrCellsMean],3), 'FontSize',32, 'ColorbarVisible','off');

s=struct(h);
s.XAxis.TickLabelRotation = 0;
%s.TitleHandle.FontWeight = 'normal';
s.Axes.XAxisLocation = 'top';
s.XAxis.TickLabelInterpreter = 'latex';
h.YDisplayLabels=repmat({''},length(s.YDisplayLabels),1);

saveas(gcf, [savePrefix 'MARD-MeanCellsR-heatmap' '.png']);