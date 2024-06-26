function genExAnalysis(cellTable,MNodesCell, options)

isRndTest=options.isRndTest;
isPlotPDF=options.isPlotPDF;
isAddNoise=options.isAddNoise;
pvalAllHMapMin=abs(log10(options.pvalAllHMapMin));
outputFolderName=options.outputFolderName;

isTwoSided=options.isTwoSided;
W=options.W;


NIDs=cellTable.NID;

cTypes=cellTable.nodeType;

minPvalEval=options.pvalAllHMapMin;



miNumCell=5;
minPVolcan=options.pvalAllHMapMin;
minFoldChange=0;
pvalAllTheoMin=3;
tablePvalQMin=0.8;

geneAnnotes=cellTable.Properties.VariableNames(options.gStart:end);
geneExpression=cellTable(:, options.gStart:end);
if ~isempty(geneExpression)
geneExpression=table2array(geneExpression);
if isAddNoise
    geneExpression=geneExpression+1e-5*rand(size(geneExpression)).*(geneExpression>0);
end

allNodeIdx=(1:length(cTypes));


foldChangeRandC={};
pvalueMedRandC={};

iResult=1;
sResults=0;



pvalueMedC=cell(1000, 1);
pvalueMedTheoC=cell(1000, 1);
foldChangeC=cell(1000, 1);
addressC=cell(1000, 1);

medMotC=cell(1000, 1);

medtypeC=cell(1000, 1);

numMotCellsC=cell(1000, 1);

geneSet=cell(1000, 1);
numExtMotifs=length(MNodesCell);

for iMt=1:numExtMotifs
    fprintf('Evaluating gene expression for motif number:  %d\n', iMt);
    motifPaths=MNodesCell{iMt};
    motifPathsType=cTypes(motifPaths);
    motifPathsTypeRv=motifPathsType(:, end:-1:1);

    if all(motifPathsType(:)==motifPathsTypeRv(:))

        postLoop=(1:ceil(W/2));
    else
        postLoop=(1:W);
    end

    for iMTPosti=postLoop
        motifPathsFirst=motifPaths(:, iMTPosti);
        motifPathsFirst=motifPathsFirst(:);
        motifPathsFirst=unique(motifPathsFirst);
        AMotifPathsTypes=cTypes(motifPathsFirst);
        [AMotifPathsTypesU, ~, jPU]=unique(AMotifPathsTypes);
        cntPT=accumarray(jPU, 1);
        [cntPTS, ist]=sort(cntPT);
        AMotifPathsTypesU=AMotifPathsTypesU(ist);
        AMotifPathsTypesU=AMotifPathsTypesU(cntPTS>miNumCell);
        cntPTS=cntPTS(cntPTS>miNumCell);
        for iMType=1:length(AMotifPathsTypesU)

            geneExpressionL=geneExpression;


            AMotifPathsi=motifPathsFirst(AMotifPathsTypes==AMotifPathsTypesU(iMType));

            typeNodeIdx=allNodeIdx(cTypes==AMotifPathsTypesU(iMType));

            motifInType=ismember(typeNodeIdx, AMotifPathsi);

            if isRndTest
                MAnimal_IDs=NIDs(AMotifPathsi(:, 1));
                TAnimal_IDs=NIDs(typeNodeIdx);

                motifInTypeRand=false(size(motifInType));

                MAnimal_IDsU=unique(MAnimal_IDs);
                for iMAI=1:length(MAnimal_IDsU)
                    flagTi=TAnimal_IDs==MAnimal_IDsU(iMAI);
                    motifInTypei=motifInType(flagTi);
                    motifInTypeRand(flagTi)=motifInTypei(randperm(length(motifInTypei)));
                end
                motifInType=motifInTypeRand;
            end

            AMotifPathsAllFlagi=(cTypes==AMotifPathsTypesU(iMType));
            geneExAlli=geneExpressionL(AMotifPathsAllFlagi, :);
            geneExi=geneExAlli(motifInType, :);
            geneNExi=geneExAlli;
            numMotCells=[sum(motifInType),sum(~motifInType)];
            medMot=median(geneExi);
            medMot=medMot(:);
            medtypeMot=median(geneNExi);
            medtypeMot=medtypeMot(:);

            [pvalueMedTheo, meDeltaMot]= computeGeneDProfile(geneExAlli, motifInType, isTwoSided);
            EvalFlag= pvalueMedTheo<minPvalEval;
            if isPlotPDF
                if any(EvalFlag)
                    tPValH=pvalueMedTheo(EvalFlag);
                    geneEvalAnnotes=geneAnnotes(EvalFlag);
                    geneExAllEval=geneExAlli(:, EvalFlag);
                    [pdfMDCell, edgeMDCell]= computeGeneDCDF(geneExAllEval, motifInType);
                    numPlots=sum(EvalFlag);
                    meDeltaMotPlot=meDeltaMot(EvalFlag);
                    motAddr=[iMt, iMTPosti, AMotifPathsTypesU(iMType),  cntPTS(iMType), length(motifInType)];

                    for iplot=1:numPlots
                        subplot(numPlots, 1,iplot)
                        pdfPlot=max(pdfMDCell{iplot}, 0);
                        edgePlot=edgeMDCell{iplot};
                        idxRight0=find(cumsum(pdfPlot)>0.999, 1);
                        idxRight1=find(edgePlot<=meDeltaMotPlot(iplot), 1, 'last');
                        idxRight0=max(idxRight0, min(idxRight1+1, length(edgePlot)));
                        tmp=find(cumsum(pdfPlot(end:-1:1))>1-eps, 1);
                        idxLeft0=length(pdfPlot);
                        if ~isempty(tmp)
                            idxLeft0=idxLeft0-tmp+1;
                        end
                        idxLeft1=length(edgePlot)-find(edgePlot(end:-1:1)<=meDeltaMotPlot(iplot), 1)+1;
                        idxLeft0=min(idxLeft0, idxLeft1);
                        if idxLeft0<5
                            idxLeft0=1;
                        end

                        pdfPlot=pdfPlot(idxLeft0:idxRight0);
                        edgePlot=edgePlot(idxLeft0:idxRight0);
                        figure('visible','off');

                        scatter(edgePlot, pdfPlot,50, 'filled')
                        hold on
                        plot(edgePlot, pdfPlot, '.-.', 'LineWidth',2, 'Color','b')
                        xlim([min(edgePlot(1), meDeltaMotPlot(iplot)-1), max(edgePlot(end), meDeltaMotPlot(iplot)+1)])
                        xline(meDeltaMotPlot(iplot),'-', {'Motif','Delta', sprintf('%2.4f',meDeltaMotPlot(iplot))}, 'color', 'r', 'linewidth', 2)
                        if iplot==1
                            titleTex=sprintf('theo pval=%2.2f with motif addr %d,%d,%d(%d/%d), gene=%s\n', abs(log(tPValH(iplot))), motAddr, string(geneEvalAnnotes{iplot}));
                        else
                            titleTex=sprintf('theo pval=%2.2f, gene=%s\n', abs(log(tPValH(iplot))), string(geneEvalAnnotes{iplot}));
                        end
                        title(titleTex);
                        grid on
                    end
                    figName=sprintf('MTAddr-%d-%d-%d',motAddr(1:3));
                    if isRndTest
                        figName=strcat('Rand', figName);
                    end
                    if ~isTwoSided
                        figName=strcat('OS', figName);
                    end

                    figName=strcat(outputFolderName, figName, '.jpg');

                    xlabel('delta median')
                    ylabel('probability')
                    saveas(gcf, figName)

                end
            end
            pvalueMedo=pvalueMedTheo;


            foldChange=meDeltaMot(:);



            foldChangeN=abs(foldChange)/max(abs(foldChange));

            siGenes=(pvalueMedo<1) & (foldChangeN>=minFoldChange);

            pvalueMed=pvalueMedo(siGenes);


            pvalueMedC{iResult}=pvalueMed;
            pvalueMedTheoC{iResult}=pvalueMedTheo(siGenes);
            foldChangeC{iResult}=foldChange(siGenes);
            addressC{iResult}=[iMt, iMTPosti, AMotifPathsTypesU(iMType)];

            medMotC{iResult}=medMot(siGenes);

            medtypeC{iResult}=medtypeMot(siGenes);

            numMotCellsC{iResult}=numMotCells;



            geneset0=(1:size(geneExpression, 2)).';

            geneSet{iResult}=geneset0(siGenes);
            iResult=iResult+1;
            sResults=sResults+length(pvalueMed);

        end
    end
end

% save Results

fprintf('Output results are being saved at %s\n', outputFolderName);


geneSet=geneSet(1:iResult-1);
pvalueMedC=pvalueMedC(1:iResult-1);
pvalueMedTheoC=pvalueMedTheoC(1:iResult-1);
foldChangeC=foldChangeC(1:iResult-1);
addressC=addressC(1:iResult-1);
medMotC=medMotC(1:iResult-1);
medtypeC=medtypeC(1:iResult-1);
numMotCellsC=numMotCellsC(1:iResult-1);



%% plot volcano


addressAll0=vertcat(addressC{:});



% heatmap for all cases, ordered with address
lengthV=zeros(length(pvalueMedC),1);
heatMapVal0=zeros(length(pvalueMedC), length(geneAnnotes));
heatMapValMD0=zeros(length(pvalueMedC), length(geneAnnotes));

for iPV=1:length(pvalueMedC)
    lengthV(iPV)=length(pvalueMedC{iPV});
    heatMapValMD0(iPV, geneSet{iPV})=foldChangeC{iPV};

    heatMapVal0(iPV, geneSet{iPV})=abs(log(pvalueMedC{iPV}));
end




heatmapYLabel=addressAll0;
[heatmapYLabel, iAddrSort]=sortrows(heatmapYLabel);
addressAll0O=addressAll0;
addressAll0=addressAll0(iAddrSort, :);
[heatmapYLabelU, iU]=unique(heatmapYLabel(:, 1));
heatmapYLabel=strcat(num2str(heatmapYLabel(:, 1)),'-', num2str(heatmapYLabel(:, 2)),'-', num2str(heatmapYLabel(:, 3)));

heatmapYLabelStr=repelem({''},length(heatmapYLabel));

iUL=iU+[iU(2:end);length(heatmapYLabel)];
iUL=round(iUL/2);

for iiU=1:length(iU)

    heatmapYLabelStr{iUL(iiU)}=num2str(heatmapYLabelU(iiU));
end

heatMapVal=heatMapVal0(iAddrSort, :);
heatMapValMD=heatMapValMD0(iAddrSort, :);



%%%%  heatmap

figure('WindowState','maximized', 'visible','off')
h=heatmap(heatMapVal,'XLabel' , 'gene','YLabel' , 'motif number','Colormap',cool,'CellLabelColor','none');



h.GridVisible = 'on';
c=gray;
c = flipud(c);
colormap(c);
Ax = gca;
Ax.XDisplayLabels = nan(length(Ax.XDisplayData),1);
Ax.YDisplayLabels = nan(length(Ax.YDisplayData),1);
% set(gca,'ColorScaling','log') % Log scale
% axs = struct(gca);
% cb = axs.Colorbar;
% 
% cb.Ticks = unique(floor(linspace(min(heatMapVal(:)), max(heatMapVal(:)), 8))) ;
% 
% 
% Ax.XDisplayLabels=geneAnnotes;

% set(struct(h).NodeChildren(3), 'XTickLabelRotation', 0); % put instead of the last example line
Ax.YDisplayLabels=heatmapYLabel;

colorbar
title('log pvalue');
figName=strcat(outputFolderName, 'allHeatmap.jpg');
saveas(gcf, figName)

%%%%%%%%%%%%%%%%%%


% heatmap for selected ones

geneFlags=any(heatMapVal>pvalAllHMapMin);
postFalgs=any(heatMapVal>pvalAllHMapMin, 2);
motifNum=addressAll0(:, 1);
numSigs=accumarray(motifNum, postFalgs);
postFalgs=numSigs>0;
postFalgs=postFalgs(motifNum);
addressAllSel=addressAll0(postFalgs, :);
heatMapValSel=heatMapVal(postFalgs, :);
heatMapValMDSel=heatMapValMD(postFalgs, :);

heatMapValSel=heatMapValSel(:, geneFlags);
heatMapValMDSel=heatMapValMDSel(:, geneFlags);

geneAnnotesSell=geneAnnotes(geneFlags);
heatmapYLabelSel=heatmapYLabel(postFalgs,:);
figure('WindowState','maximized', 'visible','off')
h=heatmap(heatMapValSel,'FontSize',12,'Colormap',cool,'CellLabelColor','none');
c=gray;
c = flipud(c);
colormap(c);
% hHeatmap = struct(h).Heatmap;
% hHeatmap.GridLineStyle = ':';

Ax = gca;
Ax.XDisplayLabels=geneAnnotesSell;
Ax.YDisplayLabels=heatmapYLabelSel;
% set(gca,'ColorScaling','log')
% axs = struct(gca);
% cb = axs.Colorbar;
% 
% cb.Ticks = unique(floor(linspace(min(heatMapVal(:)), max(heatMapVal(:)), 8))) ;


% figName=strcat(outputFolderName, 'selectedHeatmap.jpg');
% saveas(gcf, figName)

% heatmap for selected ones, ordered with type

heatmapYLabel=addressAllSel;
if ~isempty(heatmapYLabel)
    [~, iAddrSort]=sortrows([heatmapYLabel(:, end),heatmapYLabel(:, 1:end-1)] );
    heatmapYLabel=heatmapYLabel(iAddrSort, :);

    heatmapYLabel=strcat(num2str(heatmapYLabel(:, 1)),'-', num2str(heatmapYLabel(:, 2)),'-', num2str(heatmapYLabel(:, 3)));
    heatMapValSel=heatMapValSel(iAddrSort, :);
    heatMapValMDSel=heatMapValMDSel(iAddrSort, :);

else
    heatmapYLabel=heatmapYLabelSel;
end



heatmapYLabelSel=heatmapYLabel;
figure('visible','off')
h=heatmap(heatMapValSel,'FontSize',12,'Colormap',cool,'CellLabelColor','white', 'CellLabelFormat', '%.0f');
c=gray;
c = flipud(c);
colormap(c);
% hHeatmap = struct(h).Heatmap;
% hHeatmap.GridLineStyle = ':';

Ax = gca;
Ax.XDisplayLabels=geneAnnotesSell;
Ax.YDisplayLabels=heatmapYLabelSel;
% set(gca,'ColorScaling','log')
% axs = struct(gca);
% cb = axs.Colorbar;
% 
% cb.Ticks = unique(floor(linspace(min(heatMapVal(:)), max(heatMapVal(:)), 8)));
title('log pvalue');

figName=strcat(outputFolderName, 'selectedSortedHeatmap.jpg');
saveas(gcf, figName)

figure('visible','off')
h=heatmap(heatMapValMDSel,'FontSize',12,'Colormap',cool,'CellLabelColor','white', 'CellLabelFormat', '%1.1f');
c=turbo;
colormap(c);
% hHeatmap = struct(h).Heatmap;
% hHeatmap.GridLineStyle = ':';

Ax = gca;
Ax.XDisplayLabels=geneAnnotesSell;
Ax.YDisplayLabels=heatmapYLabelSel;

title('delta median');

figName=strcat(outputFolderName, 'selectedDMedianHeatmap.jpg');
saveas(gcf, figName)



pvalAllHist=vertcat(pvalueMedC{:});
pvalAllRandHist=vertcat(pvalueMedRandC{:});

pvalAllHist=abs(log(pvalAllHist));
pvalAllRandHist=abs(log(pvalAllRandHist));

isInfPval=isinf(pvalAllHist);
pvalAllHist(isInfPval)=ceil(max(pvalAllHist(~isInfPval))+1);

% [pvalAllHistSort, iHSort]=sort(pvalAllHist);


pvalAllTheo=vertcat(pvalueMedTheoC{:});

pvalAllTheo=abs(log(pvalAllTheo));


isInfPval=isinf(pvalAllTheo);
pvalAllTheo(isInfPval)=ceil(max(pvalAllTheo(~isInfPval))+1);



% pvalAllFitSort=pvalAllTheo(iHSort);


pvalAllFlag=pvalAllTheo>pvalAllTheoMin;

geneSetAll=vertcat(geneSet{:});

geneSetAllS=geneSetAll(pvalAllFlag);
pvalAllHistS=pvalAllHist(pvalAllFlag);


pvalAllFitS=pvalAllTheo(pvalAllFlag);

foldChangeAll=vertcat(foldChangeC{:});

foldChangeRandAll=vertcat(foldChangeRandC{:});



figure('visible','off')
scatter(foldChangeAll, pvalAllHist, 'filled')
hold on
scatter(foldChangeRandAll, pvalAllRandHist, 'filled')

grid minor
legend('motif gene expression pvalues')
xlabel('delta median')
ylabel('absolute log pvalue')

    figName='volcano';
    figName=strcat(outputFolderName, figName, '.jpg');

    saveas(gcf,figName);



textingFlag= pvalAllHistS>=quantile(pvalAllFitS, tablePvalQMin);
foldChangeAllS=foldChangeAll(pvalAllFlag);
addressAll=repelem(addressAll0O, lengthV,1);
numMotAllV=vertcat(numMotCellsC{:});
numMotAllV=repelem(numMotAllV, lengthV,1);
addressAllS=addressAll(pvalAllFlag, :);


geneTextD=geneSetAllS(textingFlag);
geneText=geneAnnotes(geneTextD);

geneName=geneText(:);

motifAddr=addressAllS(textingFlag, :);
motifNumber=motifAddr(:, 1);
motifPosition=motifAddr(:,2);
cellType=motifAddr(:,3);
pvalue=pvalAllFitS(textingFlag);
deltaMedian=foldChangeAllS(textingFlag);
medMotVec=vertcat(medMotC{:});
medtypeVec=vertcat(medtypeC{:});
medMotVecS=medMotVec(pvalAllFlag);
medtypeVecS=medtypeVec(pvalAllFlag);
motifMedian=medMotVecS(textingFlag);
typeMedian=medtypeVecS(textingFlag);
numMotAllVS=numMotAllV(pvalAllFlag, :);
numCellsMot=numMotAllVS(textingFlag, :);

numCellMotif=numCellsMot(:, 1);
numCellType=sum(numCellsMot, 2);







asignificantGenesTable=table(geneName,motifNumber,motifPosition,cellType, pvalue, deltaMedian, motifMedian, typeMedian, numCellMotif, numCellType);

[~, ist]=sort(pvalue, 'descend');

asignificantGenesTable=asignificantGenesTable(ist, :);



    csvFileName='geaTable';
    csvFileName=strcat(outputFolderName, csvFileName, '.csv');

    writetable(asignificantGenesTable,csvFileName)
else
    printf('There is no gene expression data!\n')

end




