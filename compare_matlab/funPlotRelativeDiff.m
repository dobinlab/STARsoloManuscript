function funPlotRelativeDiff(RD, vLine1, vColor1, lineW)

nCases=length(RD);

%% number of cell/gene vs relative devaition to ref-tool
% myFigure(201); axes;
myFigure(301); axes;    

%leg1={};
for ic1=1:nCases

    fprintf(1, '%i ', ic1);
    %leg1={leg1{:} caseNames{ic1}};      

    ma1=RD{ic1};

%     figure(201);
%     line((1:length(ma1))/length(ma1),sort(abs(ma1)), 'LineWidth',4, 'Color', vColor1{ic1}, 'LineStyle', vLine1{ic1});
%     xlabel('Proportion of expressed genes/cells')
%     ylabel('Absolute Relative Deviation')
%     %legend(leg1,'interpreter','latex', 'Location', 'NorthWest');
%     grid on
%     box on        


    figure(301);
    %line((1:length(ma1))/length(ma1),sort(ma1), 'LineWidth',4, 'Color', vColor1{ic1}, 'LineStyle', vLine1{ic1});
    ma1p=ma1; ma1p(ma1<0)=0;
    ma1m=ma1; ma1m(ma1>0)=0;

    line((1:length(ma1))/length(ma1),sort(ma1p,'descend'), 'LineWidth',lineW, 'Color', vColor1{ic1}, 'LineStyle', vLine1{ic1});
    h2=line((1:length(ma1))/length(ma1),sort(ma1m), 'LineWidth',lineW, 'Color', vColor1{ic1}, 'LineStyle', vLine1{ic1});
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); %turn off legend for the 2nd line

    xlabel('Proportion of expressed genes/cells')
    ylabel('Relative Deviation')
    %[hh,icons,plots,txt] = legend(leg1,'interpreter','latex', 'Location', 'NorthWest');
    grid on
    box on            

end    





