axpos1=get(gca,'Position');
set(gca,'Position',[axpos1(1) 0.15 axpos1(3) 0.98-0.15]);

myLegend(casesSelectNamesShort(2:end), 'SouthEast', 0.1, [-0.06 0.06], vColor1, 32);
ylabel(['Relative Deviation from ' compareReference]);
xlabel({'Proportion of expressed' 'gene/cell matrix elements'});

    
%title('Relative Deviation to true counts', 'FontWeight', 'normal', 'FontAngle',figs.titleAngle);

saveas(gcf, [savePrefix 'RD' '.png']);