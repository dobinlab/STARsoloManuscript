%% parameters used:
% caseTableFile
% casePrefix
% genome1 sample1 run1 sampleType opt_STAR opt_alevinFry
% caseLoad
% nCorr
% countThreshold
% nCells

%%
eval(caseTableFile);

if ~isempty(caseLoad)
    caseTable=caseTable(ismember(caseTable.tool, caseLoad),:); 
end

%% whitelist
if ~exist('WL','var')
    WL=funLoadCB(WLfile);
end
%%
[geneList geneListNames]=textread([casePrefix '/genomes/' genome1 '/genes.tsv'], '%s %s %*[^\n]'); %#ok<DTXTRD>
nGenes=length(geneList);

%%
% if ~exist('M','var')
%     M={};
% end
M={};
[caseTable, M]=funLoadData(M,casePrefix, caseTable, WL, geneList);

%% spliced+unslpiced for kb
vv=find(endsWith(caseTable.tool, 'total') & startsWith(caseTable.tool, 'kallisto'));
for ii=1:length(vv)
    if vv(ii)>length(M) || isempty(M{vv(ii)})
        continue;
    end
    ispl=find(strcmp(strrep(caseTable.tool{vv(ii)},'total','spliced'), caseTable.tool));
    M{vv(ii)}=M{vv(ii)}+M{ispl}; % total contained unspliced, add spliced
end

%%
casesSelect1=caseTable.tool(~cellfun(@isempty,M));

%%
[vSelect1]=ismember(caseTable.tool, casesSelect1);
caseTable1=caseTable(vSelect1,:);
nCases=size(caseTable1,1);

%%
% call nCells top cells
if nCells>0
    CBfilteredBool=funFilterTopCells(M{1},nCells); nCB=nnz(CBfilteredBool);
else
    CBfilteredBool=true(1,size(WL,1));
end

%%
Ma={[]};
[Ma, CBintersection]=funCombineMatrices(M(vSelect1), Ma, CBfilteredBool);

%%
if ~exist('corrGenesCells', 'var')
    corrGenesCells=funFullCorr(Ma, countThreshold, nCorr, caseTable1.labelSimple);
end
%
disp('Spearman')
disp([char(caseTable1.labelSimple) repmat('    ',size(corrGenesCells,1), 1) num2str(corrGenesCells(:,:,2))])

funWriteMatrixTable([savePrefix '_tables'], 'GeneCell_SpearmanR', corrGenesCells(:,:,2), caseTable1.labelSimple, caseTable1.labelSimple)

%%
if ~exist('RD', 'var')
    [RD,MARD]=funCalcRelDiff(Ma, countThreshold);
end

%%
if ~exist('corrCells', 'var')
    corrCells=funCorrBetweenCells(Ma, countThreshold);
end

%%
clear M
save([savePrefix '_Res.mat'], '-v7.3')%, 'caseTable1', 'RD', 'MARD', 'corrGenesCells', 'corrCells', 'Ma', '-v7.3');

%%
disp(['tool'                                                              ',Cell/Gene Spearman R'            ',Cell/Gene MARD'])
disp([char(caseTable1.tool(1:end)) repmat(',',nCases, 1)   num2str(corrGenesCells(1:end,1,2))   repmat(',',nCases, 1)       num2str(MARD(1:end,1))  ])

funWriteMatrixTable([savePrefix '_tables'], 'R_MARD_cellR', [ corrGenesCells(:,1,2) MARD ], caseTable1.tool, {'Cell/Gene Spearman R'   'Cell/Gene MARD'})


