function funPlotSpearmanBar(plV1, toolNames, axXstart, savePrefix, compareReference)

myFigure(51)

plN1=(length(plV1):-1:1);

barh(plN1, plV1, 0.4, 'FaceColor', [229, 209, 235]/255);
axis([axXstart 1 0.5 plN1(1)+0.5]);
set(gca, 'Position', [0.4 0.1 0.59 0.89], 'YtickLabel',toolNames, 'TickLabelInterpreter', 'none')
% xlabel(['Spearman R to ' compareReference]);
xlabel(['Spearman R'])
set(gca, 'XGrid','on');
%t1=title('Genes/cells Spearman correlation', 'HorizontalAlignment', 'right', 'Position', [1 plN1(1)+0.55 0], 'FontAngle','italic' );


plVround = round(plV1,3); % 3 digits

for ii=1:length(plVround)
    if plVround(ii)==1
        plVround(ii)=round(plV1(ii),ceil(-log10(1-plV1(1)))+1);
    end
    plText{ii}=num2str(plVround(ii),'%0.15g');
end

text(plV1, plN1, plText, 'HorizontalAlignment','right','VerticalAlignment','middle', 'Color', [0 0 0]);
saveas(gcf, [savePrefix 'SpearmanBar' '.png']);

