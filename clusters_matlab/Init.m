%%
set(0,'defaultAxesPosition',[0.15 0.12 0.8 0.8],'defaultAxesFontUnits','points','defaultAxesFontSize',28);
set(0, 'defaultFigurePaperPositionMode', 'auto','defaultFigurePosition',[50 50 1200 1200], 'defaultFigureColor', 'White');
addpath /scratch/dobin/Software/MatlabToolbox/
addpath /scratch/dobin/Software/MyMatlabToolbox/
% addpath /sonas-hs/gingeras/nlsas_norepl/user/dobin/Analysis/Common/
% addpath \\vega\dobin/data/scratch/Software/MatlabToolbox/
% addpath \\vega\dobin/data/scratch/Software/MyMatlabToolbox/

dbstop if error
% vColor={'Red','Green','Blue','Cyan','Magenta'};
% for ii=1:5
%     vColor=[vColor vColor];
% end

%%
global vLine vColor vMarker

vLine={'-','--','-.',':'};
for ii=1:5
    vLine=[vLine vLine];
end


vColor=[         
         0.514235294117647         0.769529411764706         0.725764705882353
                      0.93         0.864352941176471         0.404823529411765
         0.692941176470588         0.678352941176471         0.795058823529412
         0.915411764705882         0.466823529411765         0.415764705882353
         0.466823529411765         0.645529411764706         0.769529411764706
         0.922705882352941         0.656470588235294         0.357411764705882
         0.652823529411765         0.809647058823529         0.382941176470588
         0.685647058823529         0.466823529411765         0.689294117647059
         0.791411764705882         0.791411764705882         0.791411764705882
                     0.744         0.857058823529412         0.718470588235294
         0.919058823529412         0.747647058823529         0.835176470588235
                      0.93                      0.93         0.652823529411765
      ];
                  
set(0, 'defaultAxesColorOrder', vColor);

vMarker={'s','o','*','d','^','v','>','<'};
%%
