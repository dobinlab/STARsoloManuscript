function myLegend(legText, legLoc, shiftEnd, shiftLegPos, vColor1, legFontSize)

[hh,icons,plots,txt] = legend(legText,'interpreter','none', 'Location', legLoc, 'FontSize', legFontSize, 'FontWeight','normal');

nLeg=length(icons)/3;

for ileg=1 : nLeg
    icons(ileg).Position=icons(ileg).Position + [shiftEnd 0 0]; %text
    icons(ileg).Color=vColor1{ileg};
    icons(nLeg+(ileg-1)*2+1).XData=icons(nLeg+(ileg-1)*2+1).XData + [0 shiftEnd];%line
end
hh.Position=hh.Position+[shiftLegPos(1)   0   shiftLegPos(2)   0];