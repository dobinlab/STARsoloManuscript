run('../compare_matlab/Init.m')
addpath('../compare_matlab/');

%%
toolsID={'CR', 'Ssp', 'Sfu', 'Adf', 'Ask', 'kb'};
toolNames={'CellRanger', 'STAR_sparseSA', 'STAR_fullSA', 'Alevin_full-decoy', 'Alevin_sketch', 'Kallisto'};
toolNamesTex={'CellRanger', 'STAR\_sparseSA', 'STAR\_fullSA', 'Alevin\_full-decoy', 'Alevin\_sketch', 'Kallisto'};

dirClusters = '../clusters_scanpy/write_pbmc5k/CRclusters9/';
dirFigs = './figures_pbmc5k/CRclusters9/';
mkdir(dirFigs);


clusterNames = textread([dirClusters 'clusterTypes.txt'], '%s', 'whitespace', '\t');
nCl = length(clusterNames);

nFields = 3;
iScore=1; iPadj=2; iLogfc=3;

clear genes
[genesList] = textread([dirClusters 'genesListFull.txt'], '%s');
genesList=sort(genesList);
nGenes=length(genesList);
nTools=length(toolsID);

padjThr = 0.01;

%%
mDE=nan(nGenes, nFields, nCl, nTools);
for it = 1:nTools
    tb1 = readtable([dirClusters  toolsID{it} '.rgg_wilcoxon.csv']);
    
    for ic = 1:nCl
        g1 = table2array(tb1(:, (ic-1)*(nFields+1) + 1));
        [~, g1i] = ismember(g1, genesList);
        
        mDE(g1i,:,ic,it) = table2array(tb1(:, (ic-1)*(nFields+1) + (2:nFields+1)));
        
    end
    
end

%% figure: compare ranks
figSelect.plot = true;
figSelect.logfc = true;

figSelect.ranks = false;
figSelect.rankRD = false;
figSelect.rankJI = false; % good plot but noisy
figSelect.score = false;
figSelect.scoreDiff = false;
figSelect.scoreHist = false;
figSelect.scoreLog = false;

pearsonR_logfc=-ones(nTools, nTools, nCl);
signPosJI=zeros(nTools, nTools, nCl);
rankJI=zeros(nTools,nTools,nCl,2000);
rankUmI=zeros(nTools,nTools,nCl,2000);

padjThr = 0.01;
clear s1 s2 vv v1 v2

for ic = 1:nCl
    for it1 = 1%:length(toolsID)
        for it2 = it1+1:length(toolsID)
            s1 = mDE(:,:,ic,it1);
            s2 = mDE(:,:,ic,it2);
            
            vv = ~isnan(s1(:,1)) | ~isnan(s2(:,1));
            maxRank=nnz(vv);
            s1 = s1(vv,:);
            s2 = s2(vv,:);
            
            v1=isnan(s1(:,1));
            v2=isnan(s2(:,1));
            
            s1(v1,iPadj)=1; % p=1 for absent genes
            s2(v2,iPadj)=1;
            
            g1=s1(:,iPadj)<padjThr & s1(:,iScore)>0; % significant and positive
            g2=s2(:,iPadj)<padjThr & s2(:,iScore)>0;
            signPosJI(it1,it2,ic)=funJaccardIndex(find(g1),find(g2));
            signPosJI(it2,it1,ic)=signPosJI(it1,it2,ic);
            
%             
%             s1(v1,iLogfc)=0; % logFC=0 for absent genes
%             s2(v2,iLogfc)=0;                  
            
            %%
            if figSelect.ranks 
                s1 = mDE(:,:,ic,it1);
                s2 = mDE(:,:,ic,it2);

                vv = ~isnan(s1(:,1)) | ~isnan(s2(:,1));
                maxRank=nnz(vv);
                s1 = s1(vv,:);
                s2 = s2(vv,:);

                v1=isnan(s1(:,1));
                v2=isnan(s2(:,1));

                s1(v1,iPadj)=1; % p=1 for absent genes
                s2(v2,iPadj)=1;      

                s1(v1,iLogfc)=0; % logFC=0 for absent genes
                s2(v2,iLogfc)=0;     
                
                s1(v1,iScore)=-1e6; % low score for absent genes
                s2(v2,iScore)=-1e6;
                rank1 = tiedrank(-s1(:,iScore)); % rank descends by score
                rank2 = tiedrank(-s2(:,iScore));
                
                funScatterDots(1000+ic*100+it1*10+it2, rank1, rank2, s1(:,iPadj)<padjThr, s2(:,iPadj)<padjThr);
                
                xlabel(['Gene Rank: '   toolNames{it1}]);
                ylabel(['Gene Rank: '   toolNames{it2}]);
                
                axis([1 maxRank 1 maxRank]);                
                %set(gca,'Ytick',(0.5:0.5:2)*1e4)
                %set(gca, 'Xscale','log', 'Yscale','log');
            end

            %%
            if figSelect.ranks 
                s1 = mDE(:,:,ic,it1);
                s2 = mDE(:,:,ic,it2);

                vv = ~isnan(s1(:,1)) | ~isnan(s2(:,1));
                maxRank=nnz(vv);
                s1 = s1(vv,:);
                s2 = s2(vv,:);

                v1=isnan(s1(:,1));
                v2=isnan(s2(:,1));

                s1(v1,iPadj)=1; % p=1 for absent genes
                s2(v2,iPadj)=1;      

                s1(v1,iLogfc)=0; % logFC=0 for absent genes
                s2(v2,iLogfc)=0;     
                
                s1(v1,iScore)=-1e6; % low score for absent genes
                s2(v2,iScore)=-1e6;
                rank1 = tiedrank(-s1(:,iScore)); % rank descends by score
                rank2 = tiedrank(-s2(:,iScore));
                
                funScatterDots(1000+ic*100+it1*10+it2, rank1, rank2, s1(:,iPadj)<padjThr, s2(:,iPadj)<padjThr);
                
                xlabel(['Gene Rank: '   toolNames{it1}]);
                ylabel(['Gene Rank: '   toolNames{it2}]);
                
                axis([1 maxRank 1 maxRank]);                
                %set(gca,'Ytick',(0.5:0.5:2)*1e4)
                %set(gca, 'Xscale','log', 'Yscale','log');
            end            
            
            %%
            if figSelect.rankJI
                s1 = mDE(:,:,ic,it1);
                s2 = mDE(:,:,ic,it2);

                vv = ~isnan(s1(:,1)) | ~isnan(s2(:,1));
                maxRank=nnz(vv);
                s1 = s1(vv,:);
                s2 = s2(vv,:);

                v1=isnan(s1(:,1));
                v2=isnan(s2(:,1));

                s1(v1,iPadj)=1; % p=1 for absent genes
                s2(v2,iPadj)=1;      

                s1(v1,iLogfc)=0; % logFC=0 for absent genes
                s2(v2,iLogfc)=0;                     
                
                s1(v1,iScore)=-1e6; % low score for absent genes
                s2(v2,iScore)=-1e6;                
                rank1 = tiedrank(-s1(:,iScore)); % rank descends by score
                rank2 = tiedrank(-s2(:,iScore));
                
                [~,gr1] = sort(rank1);
                [~,gr2] = sort(rank2);
                
                for irr=1:size(rankJI,4)
                    [rankJI(it1,it2,ic,irr), rankUmI(it1,it2,ic,irr)]=funJaccardIndex(gr1(1:irr),gr2(1:irr));
                end
                
                %set(gca, 'Yscale','log');
            end
            %%
            if figSelect.score
                s1 = mDE(:,:,ic,it1);
                s2 = mDE(:,:,ic,it2);

                vv = ~isnan(s1(:,1)) | ~isnan(s2(:,1));
                maxRank=nnz(vv);
                s1 = s1(vv,:);
                s2 = s2(vv,:);

                v1=isnan(s1(:,1));
                v2=isnan(s2(:,1));

                s1(v1,iPadj)=1; % p=1 for absent genes
                s2(v2,iPadj)=1;      

                s1(v1,iLogfc)=0; % logFC=0 for absent genes
                s2(v2,iLogfc)=0;                     
                
                s1(v1,iScore)=0; % low score for absent genes
                s2(v2,iScore)=0;           
                
                s1(v1,iScore)=0; % 0 score for absent genes
                s2(v2,iScore)=0;
                
                funScatterDots(2000+ic*100+it1*10+it2, s1(:,iScore), s2(:,iScore), s1(:,iPadj)<padjThr, s2(:,iPadj)<padjThr);

                xlabel(['Gene Score: '   toolNames{it1}]);
                ylabel(['Gene Score: '   toolNames{it2}]);
                %set(gca, 'Xscale','log', 'Yscale','log');
                %axis([1 maxRank 1 maxRank]);
                %set(gca,'Ytick',(0.5:0.5:2)*1e4)            
            end
            
            %%
            if figSelect.scoreLog
                s1 = mDE(:,:,ic,it1);
                s2 = mDE(:,:,ic,it2);

                vv = ~isnan(s1(:,1)) | ~isnan(s2(:,1));
                maxRank=nnz(vv);
                s1 = s1(vv,:);
                s2 = s2(vv,:);

                v1=isnan(s1(:,1));
                v2=isnan(s2(:,1));

                s1(v1,iPadj)=1; % p=1 for absent genes
                s2(v2,iPadj)=1;      

                s1(v1,iLogfc)=0; % logFC=0 for absent genes
                s2(v2,iLogfc)=0;                     
                
                minScore = -64;
                maxScore = -minScore;
                s1(v1,iScore)=minScore; % 0 score for absent genes
                s2(v2,iScore)=minScore;
                
                funScatterDots(2000+ic*100+it1*10+it2, symLog(s1(:,iScore)), symLog(s2(:,iScore)), s1(:,iPadj)<padjThr, s2(:,iPadj)<padjThr, toolNames{it1}, toolNames{it2});
                
                axis(symLog([minScore maxScore minScore maxScore]));
                
                %vtick = [-64 -32 -16 -8 -4 -2 -1 0 1 2 4 8 16 32 64];
                vtick = [-64 -16 -4 -1 0 1 4 16 64];
                set(gca, 'Ytick',symLog(vtick), 'Xtick',symLog(vtick));
                set(gca, 'YtickLabel', vtick, 'XtickLabel',vtick);
                
                xlabel(['Gene Score: \bf{'   toolNamesTex{it1}   '}'] );
                ylabel(['Gene Score: \bf{'   toolNamesTex{it2}   '}'] );
                
                title([clusterNames{ic} ':  ' toolNamesTex{it2} ' vs ' toolNamesTex{it1}], 'FontSize', 26)
                saveas(gcf, [dirFigs 'GeneScore_cl' num2str(ic) '_' toolsID{it2} '_' toolsID{it1} '.png']);
                
                %set(gca, 'Xscale','log', 'Yscale','log');
                %axis([1 maxRank 1 maxRank]);
                %set(gca,'Ytick',(0.5:0.5:2)*1e4)            
            end            
            
            %%
            if figSelect.scoreDiff
                s1(v1,iScore) = 0; % 0 score for absent genes
                s2(v2,iScore) = 0;
                
                sdiff = s2(:,iScore)-s1(:,iScore);
                smax = max(abs([s2(:,iScore),s1(:,iScore)]),[],2);
                %smax=mean([s2(:,iScore),s1(:,iScore)],2);
                
                funScatterDots(20000+ic*100+it1*10+it2, smax, sdiff, s1(:,iPadj)<padjThr, s2(:,iPadj)<padjThr);

                xlabel(['Gene Score Max: '   toolNamesTex{it1}]);
                ylabel(['Gene Score Difference: '   toolNamesTex{it2}]);
                set(gca, 'Xscale','log')%, 'Yscale','log');
                %axis([1 maxRank 1 maxRank]);
                %set(gca,'Ytick',(0.5:0.5:2)*1e4)            
            end            
            
            %%
            if figSelect.logfc
                s1 = mDE(:,:,ic,it1);
                s2 = mDE(:,:,ic,it2);

                vv = ~isnan(s1(:,1)) | ~isnan(s2(:,1));
                maxRank=nnz(vv);
                s1 = s1(vv,:);
                s2 = s2(vv,:);

                v1=isnan(s1(:,1));
                v2=isnan(s2(:,1));

                s1(v1,iPadj)=1; % p=1 for absent genes
                s2(v2,iPadj)=1;      

                s1(v1,iLogfc)=0; % logFC=0 for absent genes
                s2(v2,iLogfc)=0;                       
                
                
                maxLog = 10;
                vvv = s1(:,iPadj)<padjThr | s2(:,iPadj)<padjThr;
                
                l1 = s1(vvv,iLogfc);
                l2 = s2(vvv,iLogfc);
                l1(l1>maxLog) = maxLog; l1(l1<-maxLog) = -maxLog;
                l2(l2>maxLog) = maxLog; l2(l2<-maxLog) = -maxLog;
                pearsonR_logfc(it1,it2,ic) = corr(l1,l2,'type','Pearson');
                pearsonR_logfc(it2,it1,ic) = pearsonR_logfc(it1,it2,ic);

                if figSelect.plot && it1==1
                
                    %vColor3 = num2cell(linspecer(3),2); %{[0 1 0], [1 0 0], [0 0 1]}
                    %vColor3 = vColor3([3 2 1]);
                    vColor3 =[0.2       0.7490    0.2
                              0.3       0.3    0.9                        
                              0.9153    0.2816    0.2878
                              ];
                    vColor3 = num2cell(vColor3,2);
                    
                    funScatterDots(3000+ic*100+it1*10+it2, l1, l2, s1(vvv,iPadj)<padjThr, s2(vvv,iPadj)<padjThr,vColor3);
                    set(gca,'FontSize',32)
                    set(gca, 'Ytick',-10:5:10, 'Xtick',-10:5:10)            

                    myLegend({'Significant in both', ['Significant in ' toolNames{it1}],  ...
                             ['Significant in ' toolNames{it2}]}, 'NorthWest', 0, [0 0], vColor3, 30);
                    title([clusterNames{ic} ':  {\bf' toolNamesTex{it2} '} vs ' toolNamesTex{it1}], 'FontSize', 26, 'FontWeight', 'normal')
                    %xlabel(['Gene log2fc: \bf{'   toolNamesTex{it1}   '}'] );
                    %ylabel(['Gene log2fc: \bf{'   toolNamesTex{it2}   '}'] );
                    xlabel(['Gene log2fc: '   toolNamesTex{it1} ] );
                    ylabel(['Gene log2fc: '   toolNamesTex{it2} ] );

                    saveas(gcf, [dirFigs 'GeneLogfc_cl' num2str(ic) '_' toolsID{it2} '_' toolsID{it1} '.png']);
                    saveas(gcf, [dirFigs 'GeneLogfc_cl' num2str(ic) '_' toolsID{it2} '_' toolsID{it1} '.eps'], 'epsc');

                    
                end
            end            
            
        end
    end
end

%% heatmap for signPosJI
myFigure(601)
set(gcf,'Position', [50 50 2000 1200]);

h=heatmap(clusterNames,  toolNames(2:end), round(squeeze(signPosJI(1,2:end,:)), 2), 'FontSize', 32, 'ColorbarVisible','off');
set(gca, 'Position', [0.26 0.26 0.73 0.69]);
title('Marker genes: Jaccard index vs CellRanger')

s=struct(h);
s.XAxis.TickLabelRotation = 45;
%h.Colormap = h.Colormap(end:-1:1,:);
s.TitleHandle.FontWeight = 'normal';
s.TitleHandle.FontSize = 36;
s.XAxis.FontSize=36;
s.YAxis.FontSize=36;

saveas(gcf, [dirFigs 'pbmc5k_cl' num2str(nCl) '_signPosJI' '.png']);

%% heatmap for pearsonR_logfc
myFigure(602)
set(gcf,'Position', [50 50 2000 1200]);

h=heatmap(clusterNames,  toolNames(2:end), round(squeeze(pearsonR_logfc(1,2:end,:)), 2), 'FontSize', 32, 'ColorbarVisible','off');
set(gca, 'Position', [0.26 0.26 0.73 0.69]);
caxis([0.7 1]);

title('Gene log2fc correlation with CellRanger ')

s=struct(h);
s.XAxis.TickLabelRotation = 45;
%h.Colormap = h.Colormap(end:-1:1,:);
s.TitleHandle.FontWeight = 'normal';
s.TitleHandle.FontSize = 36;
s.XAxis.FontSize=36;
s.YAxis.FontSize=36;

saveas(gcf, [dirFigs 'pbmc5k_cl' num2str(nCl) '_PearsonRlogfc' '.png']);

