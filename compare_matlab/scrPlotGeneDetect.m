axpos1=get(gca,'Position');
set(gca,'Position',[axpos1(1:3) 0.98-axpos1(2)]);
set(gca,'Ytick', 0:1000:10000);

myLegend(casesSelectNamesShort(2:end), 'NorthEast', 0.1, [-0.08 0.08], vColor1, 32);
xlabel(['Number of ' compareFPgenes ' genes per cell']);
ylabel('Number of cells')
%title('False positive genes per cell ', 'FontWeight', 'normal', 'FontAngle',figs.titleAngle);

saveas(gcf, [savePrefix 'GeneDetectHistFP' '.png']);

% %%
% figure(2200)
% axis('auto');xlim([0 40]);
% myLegend(casesSelectNamesShort(2:end), 'NorthEast', 0.05, [-0.03 0.03], vColor1, 32);
% xlabel('Number of False Negative genes');
% ylabel('Number of cells')
% title('Histogram of False Negative genes per cell ', 'FontWeight', 'normal');
% saveas(gcf, [savePrefix '_GeneDetectHistFN' '.png']);