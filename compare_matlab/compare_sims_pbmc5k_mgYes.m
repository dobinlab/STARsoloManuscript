Init
casePrefix='/scratch/dobin/STAR/STARsoloPreprint/maia1/';

% figTablesDir = 'FigTables/Sims-humanCR300-pbmc5k-MultiGeneYes/';
% savePrefix = [ figTablesDir ];
% load([savePrefix 'Res.mat'])%, 'caseTable1', 'RD', 'MARD', 'corrGenesCells', 'corrCells');

%% cases
genome1='human_CR_3.0.0';  sample1='pbmc_5k_sims_MultiGeneYes'; run1='/20/run3/'; sampleType='/10X/3/'; opt_STAR={}; opt_alevinFry='';
figTablesDir = 'FigTables/Sims-humanCR300-pbmc5k-MultiGeneYes/';
caseTableFile = 'caseTable_sims';
WLfile=[casePrefix 'data/whitelists/10Xv3'];

mkdir(figTablesDir)
savePrefix = [ figTablesDir ];

caseLoad = {};
nCells = 0; % 5000 for sims;  =0 for real data assuming the M{1} contains filtered data

compareReference = 'the Truth';
compareFPgenes = 'false positive';
%% parameters
% nCorr=size(Ma,3);
nCorr=0; %=0 for full correlation matrix
countThreshold = 0.5;

%% load data and do all calculations
common_LoadCalc;

%load([savePrefix 'Res.mat'])%, 'caseTable1', 'RD', 'MARD', 'corrGenesCells', 'corrCells');

%%
selectionName='bestOpts'; 
casesSelect = {'Simulated Truth' 'STARsolo fullSA mult:EM ENCODE' 'alevin-fry full-decoy cr-like-em' 'alevin-fry partial-decoy cr-like-em' 'alevin-fry sel-align cr-like-em'  'alevin-fry sketch cr-like-em' 'kallisto-bustools mult:Yes'};
casesSelectNames = casesSelect;
casesSelectNamesShort = {'Simulated Truth' 'STARsolo' 'Alevin_full-decoy' 'Alevin_partial-decoy' 'Alevin_sel-align'  'Alevin_sketch' 'Kallisto'};

[~,vSelect]=ismember(casesSelect, caseTable1.tool);

caseTableSelect=caseTable1(vSelect,:);
nCases=nnz(vSelect);
assert(nCases == length(casesSelect));
%
vSelect1 = vSelect(2:end);
vColor1=funToolsColors_Main(casesSelectNames(2:end));
vLine1=funToolsLines_Main(casesSelectNames(2:end));

%% bar plot: genes/cell SpR
corrTruth = corrGenesCells(vSelect(1), vSelect1, 2);
funPlotSpearmanBar(corrTruth, casesSelectNamesShort(end:-1:2), 0.5, savePrefix, compareReference);

%% Overall Spearman heatmap
scrPlotGenesCellsRheatmap;

%% Relative Difference
funPlotRelativeDiff(RD(vSelect1), vLine1, vColor1, 5);
%%
figure(301);
axis([0 0.41 -1 1])
scrPlotRelativeDiff;

%% Correlations between Cells Histogram
funPlotCorrBetweenCells(corrCells(:,:,vSelect1), vLine1, vColor1, 5);
%
figure(502);
axis('auto');
xlim([0.6 1]);
scrPlotCorrBetweenCells;

%% Gene detection
genesFNPmean=funPlotGeneDetect(Ma{1}, Ma(vSelect1), countThreshold, vColor1, vLine1, -25:50:2000, -2.5:5:200, 5);
funWriteMatrixTable([savePrefix '_tables'], 'FPNgenesPerCell', genesFNPmean, casesSelectNamesShort(2:end), {'False Negative genes', 'false positive genes', '% of False Negative genes', '% of false positive genes'})
%%
figure(2100)
axis('auto');
xlim([0 1500]);
scrPlotGeneDetect;

%%
scrHeatmap_MeanCellsR_MARD_meanFPgenes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Top tools selection: Relative Difference
selectionName='topToolsMult'; 
casesSelect = {'Simulated Truth' 'STARsolo fullSA mult:EM ENCODE' 'STARsolo fullSA mult:No ENCODE' ...
               'alevin-fry full-decoy cr-like-em' 'alevin-fry full-decoy cr-like'};
casesSelectNamesShort = {'Simulated Truth' 'STARsolo_mult:Yes' 'STARsolo_mult:No' ...
                         'Alevin_full-decoy_mult:Yes' 'Alevin_full-decoy_mult:No'};

[~,vSelect]=ismember(casesSelect, caseTable1.tool);

caseTableSelect=caseTable1(vSelect,:);
nCases=nnz(vSelect);
assert(nCases == length(casesSelect));
%
vSelect1 = vSelect(2:end);
%vColor1={[208,28,139]/255 [241,182,218]/255 [94,60,153]/255 [178,171,210]/255};
vColor1={[208,28,139]/255 [230,97,1]/255 [94,60,153]/255 [1,133,113]/255};
%vLine1=funToolsLines_Main(casesSelect(2:end));
vLine1={'-' '-' ':' ':'};

%%
funPlotRelativeDiff(RD(vSelect1), vLine1, vColor1, 5);
%
figure(301);
axis([0 0.07 -1 1])
axpos1=get(gca,'Position');
set(gca,'Position',[axpos1(1) 0.15 axpos1(3) 0.98-0.15]);

myLegend(casesSelectNamesShort(2:end), 'NorthEast', 0.03, [-0.03 0.03], vColor1, 32);

ylabel('Relative Deviation from the Truth');
xlabel({'Proportion of expressed' 'gene/cell matrix elements'});

% title('Relative Deviation from Truth', 'FontWeight', 'normal');

saveas(gcf, [savePrefix 'RD-' selectionName '.png']);


%% Best+Worst tools selection: MARD / mean-Spearman-cells table
selectionName='multYesVsNo'; 
casesSelect = {'Simulated Truth' 'STARsolo fullSA mult:EM ENCODE' 'STARsolo fullSA mult:No ENCODE' ...
               'alevin-fry full-decoy cr-like-em' 'alevin-fry full-decoy cr-like' 'alevin-fry sketch cr-like-em' 'alevin-fry sketch cr-like' ...
                'kallisto-bustools mult:Yes' 'kallisto-bustools mult:No'...
                };
            
casesSelectNamesShort = {'Simulated Truth' 'STARsolo_mult:Yes' 'STARsolo_mult:No' ...
                         'Alevin_full-decoy_mult:Yes' 'Alevin_full-decoy_mult:No' 'Alevin_sketch_mult:Yes' 'Alevin_sketch_mult:No' ...
                         'Kallisto_mult:Yes' 'Kallisto_mult:No'};

[vSelect, ind1]=ismember(caseTable1.tool, casesSelect);
vSelect=find(vSelect);
vSelect(nonzeros(ind1))=vSelect;
caseTableSelect=caseTable1(vSelect,:);
nCases=nnz(vSelect);
assert(nCases == length(casesSelect));
%
vSelect1 = vSelect(2:end);
vColor1=linspecer(length(vSelect1),'sequential'); vColor1 = num2cell(vColor1, 2);
vLine1=funToolsLines_Main(casesSelectNamesShort(2:end));

%corrCellsMean = squeeze(mean(corrCells(:,2,vSelect)));

MARDselect = MARD(vSelect);
label1=casesSelectNamesShort;
myFigure(601)
%h=heatmap({'MARD', 'Mean Cells R'}, label1, round([MARDselect corrCellsMean],3), 'FontSize',24, 'ColorbarVisible','off');
axes
set(gca, 'Position', [0.5 0.02 0.23 0.88]);
h=heatmap({'\bf MARD'}, label1, round([MARDselect],3), 'FontSize',28, 'ColorbarVisible','off');
axes
s=struct(h);
s.Axes.XAxisLocation = 'top';
h.Colormap = h.Colormap(end:-1:1,:);
s.XAxis.TickLabelInterpreter = 'latex';
% s.YAxis.FontSize=28;

set(gca, 'Position', [0.74 0.02 0.23 0.88]);
h=heatmap({'\bf Genes/Cells R'}, label1, round(corrGenesCells(vSelect,1,2),3), 'FontSize',28, 'ColorbarVisible','off');

s=struct(h);
s.XAxis.TickLabelRotation = 0;
%s.TitleHandle.FontWeight = 'normal';
s.Axes.XAxisLocation = 'top';
s.XAxis.TickLabelInterpreter = 'latex';
h.YDisplayLabels=repmat({''},length(s.YDisplayLabels),1);

saveas(gcf, [savePrefix 'Heatmap-MARD-GenesCellsR-' selectionName '.png']);
