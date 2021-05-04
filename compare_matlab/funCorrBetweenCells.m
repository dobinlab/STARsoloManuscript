function mCorr=funCorrBetweenCells(Ma, countThreshold, geneInd)

nCases = length(Ma);

if ~exist('geneInd','var') || isempty(geneInd)
    geneInd=1:size(Ma{1},1);
end

%% correlation within cells
mCorr=ones(size(Ma{1},2),2,nCases);
for icase=2:nCases
    mCorr(:,:,icase)=funCorrBetweenCells_FiltGenesAllCells(full(Ma{1}(geneInd,:)), full(Ma{icase}(geneInd,:)), countThreshold);
    
    fprintf(1,"%i %g %g ",icase, mean(mCorr(:,:,icase),1));
end
