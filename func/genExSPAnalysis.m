function gExResults=genExSPAnalysis(geneExpression,MNodesCell, specs)

% sparse gene expression analysis

isRndTest=specs.isRndTest;
isPlotPDF=specs.isPlotPDF;
outputFolderName=specs.outputFolderName;

isTwoSided=specs.isTwoSided;
W=specs.W;

if specs.isGEByTissue
    csTypes=[specs.cellType, specs.SID];

    [csTypesU, ~, cTypes]=unique(csTypes, 'rows');
else
    cTypes=specs.cellType;
    [csTypesU, ~, cTypes]=unique(cTypes);
end
minPvalEval=specs.pvalAllHMapMin;
miNumCell=5;
minFoldChange=0;


geneAnnotes=specs.geneAnnotes;
if ~isempty(geneExpression)

    allNodeIdx=(1:length(cTypes));

    iResult=1;
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
        % 
        % if all(specs.cellType(motifPaths(:))==specs.cellType(motifPaths(1)))
        % 
        %     postLoop=postLoop(:);
        % end
        % postLoop=postLoop(:);

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



                AMotifPathsi=motifPathsFirst(AMotifPathsTypes==AMotifPathsTypesU(iMType));
                

                if specs.motifVsAll
                    typeNodeIdx=allNodeIdx;
                    AMotifPathsAllFlagi=true(length(cTypes),1);
                else
                    typeNodeIdx=allNodeIdx(cTypes==AMotifPathsTypesU(iMType));
                    AMotifPathsAllFlagi=(cTypes==AMotifPathsTypesU(iMType));
                end

                motifInType=ismember(typeNodeIdx, AMotifPathsi);

                if isRndTest
                    % SIDs=cellTable.SID;
                    % MAIDs=SIDs(AMotifPathsi(:, 1));
                    % TAIDs=SIDs(typeNodeIdx);
                    % 
                    % motifInTypeRand=false(size(motifInType));
                    % 
                    % MAnimal_IDsU=unique(MAIDs);
                    % for iMAI=1:length(MAnimal_IDsU)
                    %     flagTi=TAIDs==MAnimal_IDsU(iMAI);
                    %     motifInTypei=motifInType(flagTi);
                    %     motifInTypeRand(flagTi)=motifInTypei(randperm(length(motifInTypei)));
                    % end
                    % motifInType=motifInTypeRand;
                    motifInType=motifInType(randperm(length(motifInType)));
                end

                
                geneExAlli=geneExpression(AMotifPathsAllFlagi, :);
                geneExi=geneExAlli(motifInType, :);
                geneNExi=geneExAlli;
                numMotCells=[sum(motifInType),sum(~motifInType)];
                medMot=median(geneExi);
                medMot=medMot(:);
                medtypeMot=median(geneNExi);
                medtypeMot=medtypeMot(:);

                [pvalueMedTheo, meDeltaMot]= computeGeneDProfile(geneExAlli, motifInType, minPvalEval);
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
                        f=figure('visible','off');
                        f.Position(4)=2*f.Position(4);



                        for iplot=1:numPlots
                            subplot(numPlots, 1,iplot)
                            
                            pdfPlot=max(pdfMDCell{iplot}, 0);
                            pdfPlot=pdfPlot/sum(pdfPlot);
                            edgePlot=edgeMDCell{iplot};
                            if length(pdfPlot)>50

                                temp=cumsum(pdfPlot);
                                temp(end)=1;



                                idxRight0=find(temp>=1-eps, 1);  
                                
                                tmp=[0; pdfPlot(:)];
                                idxLeft0=find(cumsum(tmp)==0,1, 'last')-3;
                                if idxLeft0<5 
                                    idxLeft0=1;
                                end
                  
    
                                pdfPlot=pdfPlot(idxLeft0:idxRight0);
                                edgePlot=edgePlot(idxLeft0:idxRight0);
                            end
                            


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
                        close(gcf);

                    end
                end


                foldChange=meDeltaMot(:);



                foldChangeN=abs(foldChange)/max(abs(foldChange));

                siGenes=(pvalueMedTheo<minPvalEval) & (foldChangeN>=minFoldChange);



                pvalueMedTheoC{iResult}=pvalueMedTheo(siGenes);
                foldChangeC{iResult}=foldChange(siGenes);
                addressC{iResult}=[iMt, iMTPosti(1), csTypesU(AMotifPathsTypesU(iMType), :)];

                medMotC{iResult}=medMot(siGenes);

                medtypeC{iResult}=medtypeMot(siGenes);

                numMotCellsC{iResult}=numMotCells;



                geneset0=(1:size(geneExpression, 2)).';

                geneSet{iResult}=geneset0(siGenes);
                iResult=iResult+1;

            end
        end
    end

    % save Results

    fprintf('Output results are being saved at %s\n', outputFolderName);


    geneSet=geneSet(1:iResult-1);
    pvalueMedTheoC=pvalueMedTheoC(1:iResult-1);
    foldChangeC=foldChangeC(1:iResult-1);
    addressC=addressC(1:iResult-1);
    medMotC=medMotC(1:iResult-1);
    medtypeC=medtypeC(1:iResult-1);
    numMotCellsC=numMotCellsC(1:iResult-1);


    gExResults.geneSet=geneSet;
    gExResults.pvalueMedTheoC=pvalueMedTheoC;
    gExResults.foldChangeC=foldChangeC;
    gExResults.addressC=addressC;
    gExResults.medMotC=medMotC;
    gExResults.medtypeC=medtypeC;
    gExResults.numMotCellsC=numMotCellsC;



    
else
    printf('There is no gene expression data!\n')

end




