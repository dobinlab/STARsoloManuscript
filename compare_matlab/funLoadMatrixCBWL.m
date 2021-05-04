function matOut=funLoadMatrixCBWL(prefix1, dir1, matIn, cbIn, WL, nGenes, switchRowCol, geneFile, geneList)

if ~isempty(cbIn)
    cb1=funLoadCB([prefix1 cbIn]);
    [~,cbInd]=ismember(cb1,WL,'rows');
end

if exist('geneFile','var') && ~isempty(geneFile)
    geneID=textread([prefix1 dir1 geneFile],'%s'); %#ok<DTXTRD>
    [~,geneInd]=ismember(regexprep(geneID,'\..*', ''), regexprep(geneList,'\..*', ''));
end

%%
f1=fopen([prefix1 matIn]);
m1=textscan(f1,'%f %f %f','CollectOutput',1, 'CommentStyle', '%');
fclose(f1);
m1=m1{1}(2:end,:); % skip the mat-info line

if (switchRowCol)
    m1=m1(:,[2 1 3]);
end

m1=m1(m1(:,1)<=nGenes,:);

%%
for ii=3:size(m1,2)
    if ~isempty(cbIn)
        m2=[m1(:,1) cbInd(m1(:,2)) m1(:,ii)];
		m2=m2(m2(:,2)>0,:); % in case cb1 contains barcodes not in the WL
    else
        m2=m1;
    end
	
	matOut{ii-2}=sparse(m2(:,1), m2(:,2), m2(:,3), nGenes, size(WL,1));
end

if exist('geneID','var')
    for ii=1:length(matOut)
        matOut{ii}(geneInd,:)=matOut{ii};
    end
end

if length(matOut)==1
	matOut=matOut{1};
end
%%
end
