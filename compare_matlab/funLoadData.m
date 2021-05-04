function [caseTable, M]=funLoadData(M, casePrefix, caseTable, WL, geneList)

nGenes = length(geneList);
%% load all cases. Assumes genes are always listed in the same order, i.e. have same numerical indexes
for ii=1:length(caseTable.mat)
    if exist('M') && ii<=length(M) && ~isempty(M{ii})
        continue;
    end
    if ~exist([casePrefix caseTable.mat{ii}])
        disp(['File not found: ' casePrefix caseTable.mat{ii}]);
        continue;
    end
    M{ii}=funLoadMatrixCBWL(casePrefix, caseTable.dir{ii}, caseTable.mat{ii}, caseTable.cb{ii}, WL, nGenes, caseTable.rowcol(ii), caseTable.geneFile{ii}, geneList);
    
    fprintf(1,"%i %s\n",ii,caseTable.mat{ii});
    
end

