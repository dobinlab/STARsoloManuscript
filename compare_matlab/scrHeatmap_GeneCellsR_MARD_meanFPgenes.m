%% MARD / mean-Spearman-cells table
genesCellsR=corrGenesCells(vSelect,vSelect(1),2);
corrCellsMean = squeeze(mean(corrCells(:,2,vSelect)));
MARDselect = MARD(vSelect);
FPgenes = [0; genesFNPmean(:,4)]*100;


label1=casesSelectNamesShort;

yS=0.02;
yW=0.8;

xS=0.4;
xB=0.01;
xW=(1-xS-2*xB-0.02)/3;

myFigure(601)

%%
axes
set(gca, 'Position', [xS yS xW yW]);
h=heatmap({'Genes/ Cells R'}, label1, round(genesCellsR,3), 'FontSize',32, 'ColorbarVisible','off');

s=struct(h);
%s.Axes.XAxisLocation = 'top';
%s.XAxis.TickLabelInterpreter = 'latex';
%s.XAxis.TickLabelRotation = 45;

h.XDisplayLabels={''};
annotation(gcf, 'textbox', [xS yS+yW xW 0.1], 'String', 'Genes/ Cells R', 'LineStyle', 'none', ...
            'FontSize',32, 'HorizontalAlignment', 'Center', 'FontWeight','bold');


%%
axes
set(gca, 'Position', [xS+xW+xB yS xW yW]);
h=heatmap({'MARD'}, label1, round(MARDselect,3), 'FontSize',32, 'ColorbarVisible','off');
h.Colormap = h.Colormap(end:-1:1,:);

s=struct(h);
s.XAxis.TickLabelRotation = 0;
s.Axes.XAxisLocation = 'top';
s.XAxis.FontWeight='bold';
%s.XAxis.TickLabelInterpreter = 'latex';
h.YDisplayLabels=repmat({''},length(s.YDisplayLabels),1);

%%
% arrayfun(@(x){sprintf('%0.2f',x)}, c_small);

axes
set(gca, 'Position', [xS+xW+xB+xW+xB yS xW yW]);
h=heatmap({'Mean Cells R'}, label1, round(FPgenes,1), 'FontSize',32, 'ColorbarVisible','off');
h.Colormap = h.Colormap(end:-1:1,:);
h.CellLabelFormat='%0.1f%%';

s=struct(h);
s.XAxis.TickLabelRotation = 0;
s.Axes.XAxisLocation = 'top';
s.XAxis.TickLabelInterpreter = 'latex';
h.YDisplayLabels=repmat({''},length(s.YDisplayLabels),1);

h.XDisplayLabels={''};
annotation(gcf, 'textbox', [xS+xW+xB+xW+xB yS+yW xW 0.12], 'String', 'FP genes per cell', 'LineStyle', 'none', ...
            'FontSize',32, 'HorizontalAlignment', 'Center', 'FontWeight','bold','VerticalAlignment', 'bottom');
%%
saveas(gcf, [savePrefix 'Heatmap-GenesCellsR-MARD-meanFPgenes' '.png']);