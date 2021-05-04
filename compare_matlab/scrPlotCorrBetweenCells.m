axpos1=get(gca,'Position');
set(gca,'Position',[axpos1(1:3) 0.98-axpos1(2)]);
set(gca,'Ytick', 0:1000:10000);

myLegend(casesSelectNamesShort(2:end), 'NorthWest', 0.1, [0 0.08], vColor1, 32);
xlabel(['Per-cell Spearman R']);
ylabel('Number of cells')
%title('Spearman R vs true counts in each cell', 'FontWeight', 'normal', 'FontAngle',figs.titleAngle);

saveas(gcf, [savePrefix 'CellCorrelationsSpearman' '.png']);
