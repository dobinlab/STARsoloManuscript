function funScatterDots(iFig, pl1, pl2, padj1, padj2, vColor1, tool1, tool2)

if ~exist('vColor1','var') || isempty(vColor1)
    vColor1={[0 1 0], [1 0 0], [0 0 1]};
end
    

myFigure(iFig);

vv=padj1 & padj2;
%Inf to have an empty set to keep the legend
line([pl1(vv); Inf], [pl2(vv); Inf], 'LineStyle', 'none', 'Marker', 'o', 'Color', vColor1{1}, 'MarkerFaceColor',vColor1{1});

vv=padj1 & ~padj2;
line([pl1(vv); Inf], [pl2(vv); Inf], 'LineStyle', 'none', 'Marker', 'o', 'Color', vColor1{2}, 'MarkerFaceColor',vColor1{2});            

vv=~padj1 & padj2;
line([pl1(vv); Inf], [pl2(vv); Inf], 'LineStyle', 'none', 'Marker', 'o', 'Color', vColor1{3}, 'MarkerFaceColor',vColor1{3});

vv=~padj1 & ~padj2;
if nnz(vv)>0
    line(pl1(vv), pl2(vv), 'LineStyle', 'none', 'Marker', 'o', 'Color', [0.3 0.3 0.3]);
end

box on
grid on
if nargin>=7
    tool1=strrep(tool1, '_', '\_');
    tool2=strrep(tool2, '_', '\_');
    legend({'Significant in both', ['Significant in ' tool2],  ['Significant in ' tool1], 'Not significant in both'}, 'Location', 'NorthWest')
end

