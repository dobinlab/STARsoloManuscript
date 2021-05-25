Init
casePrefix='/scratch/dobin/STAR/STARsoloPreprint/maia1/';

%%
figTablesDir = 'FigTables/Real-humanCR300-pbmc5k/';
savePrefix = [ figTablesDir ];
load([savePrefix 'Res.mat'])%, 'caseTable1', 'RD', 'MARD', 'corrGenesCells', 'corrCells');

%% cases
genome1='human_CR_3.0.0';  sample1='pbmc_5k'; run1='/20/b01/'; sampleType='/10X/3/'; opt_STAR={}; opt_alevinFry='';
figTablesDir = 'FigTables/Real-humanCR300-pbmc5k/';
caseTableFile = 'caseTable_real';
WLfile=[casePrefix 'data/whitelists/10Xv3'];

mkdir(figTablesDir)
savePrefix = [ figTablesDir ];

caseLoad = {};
nCells = 0; % 5000 for sims;  =0 for real data assuming the M{1} contains filtered data

compareReference = 'CellRanger';
compareFPgenes = 'extra';

%% parameters
% nCorr=size(Ma,3);
nCorr=0; %=0 for full correlation matrix
countThreshold = 0.5;

%% load data and do all calculations
common_LoadCalc;

%%
selectionName='bestOpts'; 
casesSelect = {'CellRanger 5.0.1' 'STAR sparseSA mult:No CR4' 'STAR fullSA mult:No CR4' 'alevin-fry full-decoy cr-like' 'alevin-fry partial-decoy cr-like' 'alevin-fry sel-align cr-like'  'alevin-fry sketch cr-like' 'kallisto-bustools mult:No'};
casesSelectNames = casesSelect;
casesSelectNamesShort = {'CellRanger' 'STARsolo_sparseSA' 'STARsolo_fullSA' 'Alevin_full-decoy' 'Alevin_partial-decoy' 'Alevin_sel-align'  'Alevin_sketch' 'Kallisto'};

[~,vSelect]=ismember(casesSelect, caseTable1.tool);
caseTableSelect=caseTable1(vSelect,:);
nCases=nnz(vSelect);
%
vSelect1 = vSelect(2:end);
vColor1=funToolsColors_Main(casesSelectNames(2:end));
vLine1=funToolsLines_Main(casesSelectNames(2:end));
vLine1{2}='--';

%% bar plot: genes/cell SpR
corrTruth = corrGenesCells(vSelect(1), vSelect1, 2);
funPlotSpearmanBar(corrTruth, casesSelectNamesShort(end:-1:2), 0.6, savePrefix, compareReference);

%% Overall Spearman heatmap
scrPlotGenesCellsRheatmap;

%% Relative Difference
funPlotRelativeDiff(RD(vSelect1), vLine1, vColor1, 5);
%%
figure(301);
axis('auto');axis([0 0.4 -1 1])
scrPlotRelativeDiff;

%% Correlations between Cells Histogram
funPlotCorrBetweenCells(corrCells(:,:,vSelect1), vLine1, vColor1, 5);
%
figure(502);
axis('auto');
xlim([0.75 1]);
scrPlotCorrBetweenCells;

%% Gene detection
genesFNPmean=funPlotGeneDetect(Ma{1}, Ma(vSelect1), countThreshold, vColor1, vLine1, -25:50:2000, -2.5:5:200, 5);
funWriteMatrixTable([savePrefix '_tables'], 'FPNgenesPerCell', genesFNPmean, casesSelectNamesShort(2:end), {'False Negative genes', 'false positive genes', '% of False Negative genes', '% of false positive genes'})
%%
figure(2100)
axis('auto');
xlim([0 1200]);
scrPlotGeneDetect;

%%
scrHeatmap_MeanCellsR_MARD_meanFPgenes;

%% Correlation plots
%funPlotScatter(Ma, countThreshold, caseTableSelect.label, caseTableSelect.tool, [figuresDir sample1]);
