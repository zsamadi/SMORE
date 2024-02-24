clc
close all
if length(findall(0))>1
    delete(findall(0));
end
clearvars

% created for noise edition of control data generation
% the output results very much depend on the holdout and training dada, so
% be careful with that

cd(fileparts(which(mfilename)));
cd('..\')

addpath(genpath(pwd))

isTestMode=false;
isNoiseTest=false;


W=4;
numExtMotifs=30;

cgOptions.isJClusterID=false;
cgOptions.gMode="delaunay"; % "delaunay" or "knn", or "epsilon"
cgOptions.numNeighs=5;




shConfig.fixedNodes=[];
bregID=1;
trNumShuffleO=5;
isSampled=false;

isDumb=0;
pathLength=W;

c=clock;
c=c(1:end-1);
cStr=string(c(1));
for ic=2:length(c)-1
    cStr=strcat(cStr, '_', string(c(ic)));
end


dataFolderName='D:\nucla\P1\matIO\data\';

outputFolderName0='D:\nucla\P1\matIO\output\hypo3D\';

outputFolderName0=strcat(outputFolderName0,'W', num2str(W));

if ~isTestMode
    outputFolderName0=strcat(outputFolderName0,'_',cStr);
end


outputFolderName0=strcat(outputFolderName0, '\');

isPlotGraph=true;
gHLoptions.is3D=true;
for shuffleMode=["kernelPath","shuffle"] % ["shuffle", "kernel", "kernelPath"]
    for numFixedTypes=9
        for ambigRatio=0%[1, 12]
            close all
            if length(findall(0))>1
                delete(findall(0));
            end
            trNumShuffle=trNumShuffleO;


            isDumb=isDumb+1;

            if isDumb<1
                continue;
            end

            samPro=0.5;






            % delete(gcp('nocreate'))
            % parpool;
            rng(1750);
            %% Hyperparameters
            % Hyperparameters for Creating Graph
            RBPDiscard=false;           % Discard RBP Cells from downward analysis.

            % Hyperparameters for Random Walk
            isEraseNodes=true;
            % numFixedTypes=14;
            fixedTypes=(0:numFixedTypes);
            fixedTypes=fixedTypes(:);


            %% Creat Graph




            plotGraph=true;

            cgOptions.bregID=bregID;

            cgOptions.ambigRatio=ambigRatio;

            cgOptions.ambigType=78;

            cgOptions.isNoiseCTS=false;
            cgOptions.fixedTypes=fixedTypes;

            cgOptions.retinaAndSection=[3, 3];
            cgOptions.hypoIDs=40;

            cgOptions.plotGraph=false;
            cgOptions.isUniformWeight=true;
            cgOptions.isRandom=false;
            cgOptions.nCTypesRand=12;
            cgOptions.nNodesRand=cgOptions.nCTypesRand*840;
            cgOptions.randCoords=false;

            cgOptions.isNBreg=true;
            cgOptions.isPBreg=true;
            cgOptions.isFemale=false;
            cgOptions.isMale=true;
            cgOptions.isNaive=true;
            cgOptions.isStimu=false;
            cgOptions.ambgSWTs=(15:84);


            outputFolderName=strcat(outputFolderName0,num2str(isDumb-1), '\');

            if ~exist(outputFolderName, 'dir')
                mkdir(outputFolderName)
            end

            logFile=strcat(outputFolderName, 'runlog.txt');
            % fileID=fopen(logFile,'w');


            fprintf('evaluating motifs with shuffling method %s and numFixedTypes %d \n', shuffleMode, numFixedTypes);



            if cgOptions.isRandom
                fignamExtnd='Random.jpeg';
                filenamExtnd='Random.mat';
                matFigExtend='Random.fig';
                textExtend='Random.txt';
                xGLLim=-30;
                xGHLim=30;
                yGLLim=-25;
                yGHLim=25;

            else
                fignamExtnd='Hypo3D.jpeg';
                filenamExtnd='Hypo3D.mat';
                matFigExtend='Hypo3D.fig';
                textExtend='Hypo3D.txt';
                xGLLim=-6000;
                xGHLim=6000;
                yGLLim=-5500;
                yGHLim=3500;
            end

            % rmeConfig.fixedTypes=fixedTypes;
            % rmeConfig.gMode=cgOptions.gMode;
            % rmeConfig.shuffleMode=shuffleMode;

            % done=writeReadmeHypo(outputFolderName,  rmeConfig);




            cgOptions.xyNoiseStd=0;

            if strcmpi(shuffleMode, 'kernel')
                cgOptions.rEpsilon=100;
                cgOptions.isKNeighs=true;
            else
                cgOptions.isKNeighs=false;
                cgOptions.rEpsilon=0;

            end

            csvName=strcat(dataFolderName, 'seurat_v2');

            [G,gStruct, cellTypesOne]=creat3DGraph(csvName,outputFolderName, cgOptions);

            cPTypes=(G.Nodes.label(:, 1)).';

            gHLoptions.ctAnnot=cellTypesOne;

            if isPlotGraph
                gHLoptions.folderName=outputFolderName;
                gHLoptions.isShuffled=false;
                plotHLightGraph(G, cPTypes, gHLoptions);
            end

            cSections=(G.Nodes.label(:, 2)).';
            cSectionsU=unique(cSections);

            shConfig.shuffleMode=shuffleMode;

            if strcmpi(shuffleMode, 'noise')
                cgOptions.xyNoiseStd=3*edgeDistStd;
                [GR,gStructR, ~,  BregmaAllUR, AnimalIDIdxR, cellTableR]=creatImGraph(T,outputFolderName,cgOptions);
                cNTypes=(GR.Nodes.label(:, 1)).';
                cSectionsR=(GR.Nodes.label(:, 2)).';
                cSectionsUR=unique(cSectionsR);
                xyCoordinatesR=GR.Nodes.Coordinates;
                [nodeListR, lengthListR, nodeSectionsR]=enumerateUSRKPaths(gStructR, pathLength, [ones(1,pathLength-1), samPro]);

                 if isPlotGraph
                    gHLoptions.folderName=outputFolderName;
                    gHLoptions.isShuffled=true;
                    
                    plotHLightGraph(GR, cNTypes, gHLoptions);
                 end

            else
                shConfig.rEpsilon=cgOptions.rEpsilon;
                shConfig.fixedTypes=fixedTypes;
                shConfig.isSectShuffle=true;
                shConfig.numClusters=20;
                shConfig.cSections=cSections;
                shConfig.fixedTypes=fixedTypes;
                shConfig.nearNeighs=G.Nodes.nearNeighs;
                shConfig.numClusters=1;

                shConfig.numShuffle=trNumShuffle;
                if strcmpi(shuffleMode, "kernelPath")
                    nodePathShuffle=enumerateUSROPaths(gStruct, 8, ones(1,8));

                    shConfig.nearNeighs=nodePathShuffle;
                end

                cNTypes=getShuffleTypes(cPTypes, shConfig);

                 if isPlotGraph
                    gHLoptions.folderName=outputFolderName;
                    gHLoptions.isShuffled=true;
                    plotHLightGraph(G, cNTypes(1:length(cPTypes)), gHLoptions);
                 end

            end


            if isNoiseTest
                randPermi=randperm(length(G.Nodes.label(:, 1)));
                G.Nodes.label(:, 1)=G.Nodes.label(randPermi, 1);
            end









            %% Type Shuffling


            cPTypesU=(1:length(cellTypesOne));


            PWMT=sum(cPTypes(:)==(cPTypesU(:)).');
            PWMT=PWMT/sum(PWMT);

            figure('WindowState' ,'maximize')
            bar(cPTypesU, PWMT)
            grid on
            xlabel('cell type')
            ylabel('frequency')
            text(cPTypesU-0.5, PWMT+0.001, cellTypesOne)
            title('Primary background')
            figname=strcat(outputFolderName,'PrimaryBg', fignamExtnd);

            saveas(gcf,figname)


            PWMTNK=sum(cNTypes(:)==cPTypesU);
            PWMTNK=PWMTNK/sum(PWMTNK);
            figure('WindowState' ,'maximize')
            bar(cPTypesU, PWMTNK)
            grid on
            xlabel('cell type')
            ylabel('frequency')
            text(cPTypesU-0.5, PWMT+0.001, cellTypesOne)
            title('kernel shuffled background')

            figname=strcat(outputFolderName,'kernelShfBg', fignamExtnd);

            saveas(gcf,figname)


            alphabet=cellTypesOne;
            numCells=length(alphabet);


            PWMT(PWMT==0)=0.0001;
            PWMTNK(PWMTNK==0)=0.0001;

            PWMRatio=PWMTNK./PWMT;

            figure

            plot(cPTypesU, PWMRatio, 'LineWidth',1.5)
            grid on
            title('Freq ratio kernel\\primari')
            figname=strcat(outputFolderName,'kernelPrimBgRatio', fignamExtnd);

            saveas(gcf,figname)

            %% Sample Paths from the Graph, Rand_ESU, or Random Walk




            if ~isSampled
                isSampled=true;
                [nodeList, lengthList, nodeSections]=enumerateUSRKPathsTest(gStruct, pathLength, [ones(1,pathLength-1),samPro]);
            end
            if strcmpi(shuffleMode, 'noise')

                [nodeListR, lengthListR, nodeSectionsR]=enumerateUSRKPaths(gStructR, pathLength, [ones(1,pathLength-1), samPro]);
            end


            %% Generate Shuffled data with k-mers presereved

            timerValue=tic;
            gOptions.hFrac=0;
            gOptions.mkvOrder=0;
            gOptions.rvp=true;

            gOptions.isSectHold=true;
            gOptions.isHExclsv=true;



            gOptions.numNodes=length(cPTypes);

            nodeSectionsGen=nodeSections;
            nodeSectionsGenu=unique(nodeSectionsGen);

            numHLD=ceil(length(nodeSectionsGenu)/4);

            randiHold=nodeSectionsGenu((1:length(nodeSectionsGenu)));

            randiHold=randiHold(1:numHLD);


            nodeSectionsGen(ismember(nodeSectionsGen,randiHold))=1;
            gOptions.numGRs=1;

            [posSeq, pHoldSeq, cPTypes0]=generateSeqs(nodeList,nodeSectionsGen,cPTypes, gOptions);
            % pHoldSeq=posSeq;


            if strcmpi(shuffleMode, 'noise')
                trNumShuffle=1;
                gOptions.numGRs=1;
                nodeSectionsGenR=nodeSectionsR;
                nodeSectionsGenRu=unique(nodeSectionsGenR);

                randiHoldR=randiHold;
                nodeSectionsGenR(ismember(nodeSectionsGenR,randiHoldR))=1;

                [negSeq, nHoldSeq, cNTypes0]=generateSeqs(nodeListR,nodeSectionsGenR,cNTypes, gOptions);
                shuffleModeGMaps= "kernel";
                GMaps=GR;
            else

                gOptions.numGRs=trNumShuffle;
                [negSeq, ~, cNTypes0]=generateSeqs(nodeList,nodeSectionsGen,cNTypes, gOptions);

                nHoldSeq=pHoldSeq;
                GMaps=G;
                shuffleModeGMaps=shuffleMode;

            end
            delete(gcp('nocreate'))




            diffMotif=false;
            indSeedMode=false;
            isUBack=false;
            isErzFixNodes=false;
            isEnrich=true;
            isOLess=false;


            [extMotif,textOut, commandText, background, cPTypesOut]=mtStreme(cPTypes=cPTypes0, cNTypes=cNTypes0,cPHTypes=cPTypes,cNHTypes=cNTypes(1:length(cPTypes)),pSeq=posSeq,nSeq=negSeq,pHSeq=pHoldSeq, nHSeq=nHoldSeq, ...
                rvp=gOptions.rvp, mkvOrder=gOptions.mkvOrder, wMin=W, wMax=W,threshold=0.005, nmotifs=numExtMotifs,alphabet=alphabet, fixedTypes=fixedTypes,nearNeighs=shConfig.nearNeighs,...
                shuffleMode=shuffleModeGMaps, isUBack=isUBack, trNumShuffle=trNumShuffle, diffMotif=diffMotif, indSeedMode=indSeedMode);

            elapsedTime=toc(timerValue);

            %
            % [extMotif,textOut, commandText, background, extMotifEval]=grStreme(cPTypes=cPTypes, cNTypes=cNTypes,cPHTypes=cHTypes,cNHTypes=cNHTypes,pSeq=posSeq,nSeq=negSeq,pHSeq=pHoldSeq, nHSeq=nHoldSeq, ...
            %     rvp=gOptions.rvp, mkvOrder=gOptions.mkvOrder, wMin=W, wMax=W,threshold=0.005, nmotifs=numExtMotifs,alphabet=alphabet, isEraseNodes=isEraseNodes);

            %%
            seqData.pSeq=posSeq;
            seqData.nSeq=negSeq;
            seqData.pHSeq=pHoldSeq;
            seqData.nHSeq=nHoldSeq;
            seqData.back=background;
            seqData.cPTypes=cPTypes;
            seqData.cNTypes=cNTypes;
            seqData.cPHTypes=cPTypes;
            seqData.cNHTypes=cNTypes;
            seqData.rvp=gOptions.rvp;
            seqData.fixedNodes=[];

            erzOptions.isErzNegative=false;

            erzOptions.fixedTypes=fixedTypes;

            erzOptions.shuffleMode=shuffleMode;
            erzOptions.G=G;
            erzOptions.isErzFixNodes=isErzFixNodes;
            erzOptions.isOLess=isOLess;
            erzOptions.trNumShuffle=trNumShuffle;
            erzOptions.cPTypesInit=cPTypes;
            erzOptions.cNTypesInit=cNTypes;
            erzOptions.isErzHOut=false;


            for iEx=1:length(extMotif)

                outMotif=extMotif{iEx};

                outMotif.PWMSE=log2(eps+outMotif.PWM)-log2(background{1});


                [~, sitesErased]=eraseMotif(outMotif,  seqData, erzOptions);

                outMotif.pNodes=sitesErased.pNodes;
                outMotif.nNodes=sitesErased.nNodes;
                outMotif.pHNodes=sitesErased.pHNodes;
                outMotif.nHNodes=sitesErased.nHNodes;

                outMotif.nsites=sum(sitesErased.pSeq)+sum(sitesErased.pHSeq);

                ppHNodes=[sitesErased.pNodes(:);sitesErased.pHNodes(:)];



                extMotif{iEx}=outMotif;

                % seqData.cPTypes(sitesErased.pNodes(:))=0;
                % seqData.cNTypes(sitesErased.nNodes(:))=0;
                % seqData.cPHTypes(sitesErased.pHNodes(:))=0;
                % seqData.cNHTypes(sitesErased.nHNodes(:))=0;



            end







            bkg.PWM0=background{1};
            bkg.mkvOrder=gOptions.mkvOrder;
            filename=strcat(outputFolderName, 'zSTREME','_Hypo3D', num2str(W), '.txt');
            writeTextOutput(filename, extMotif,bkg,  textOut,commandText, elapsedTime);


            numExtMotifs=length(extMotif);

            ProfDateNum=datetime("now");

            fprintf('shuffling method %s and numFixedTypes %d finished at : %s\n', shuffleMode, numFixedTypes, ProfDateNum);


            %% Seqlogo




            for iexMotif=1:numExtMotifs

                exMotifi=extMotif{iexMotif};

                if cgOptions.isJClusterID
                    zeroAddL=84-length(cellTypesOne);

                    aa=[PWMT,0.0001*ones(1, zeroAddL)] ;
                    aa=aa/sum(aa);
                    background{1}=aa;
                    exMotifi.background=background;
                    exMotifi.PWM=[exMotifi.PWM; zeros(zeroAddL,W)];
                else
                    exMotifi.background=background;
                end
                [~, exMotifi.cSeed]=max(exMotifi.PWM);

                [WL, hfig, cTypeChars]=seqlogoHypo3D(exMotifi);

                cSeedStr = strjoin(string(exMotifi.cSeed), '_');


                figname=strcat(outputFolderName, 'mLogo',num2str(iexMotif),'_', cSeedStr, fignamExtnd);
                saveas(hfig,figname)

            end

            pvalueVec=zeros(numExtMotifs, 1);
            scoreThrVec=zeros(numExtMotifs, 1);

            for iexMotif=1:numExtMotifs

                exMotifi=extMotif{iexMotif};
                pvalueVec(iexMotif)=exMotifi.testPvalue;
                scoreThrVec(iexMotif)=exMotifi.scoreThr;

            end

            %% pvalues

            figure
            plot(pvalueVec, '-kv', 'LineWidth',1)
            grid on

            xlabel('motif index')
            ylabel('evalue')

            figname=strcat(outputFolderName,'evalueWG', num2str(W), fignamExtnd);
            saveas(gcf,figname)

            figure
            plot(sort(scoreThrVec))
            grid on

            rmeConfig=cgOptions;

            rmeConfig.cTypeChars=cTypeChars;

            rmeConfig.cellTypes=cellTypesOne;

            rmeConfig.W=W;

            rmeConfig.randiHold=randiHold;
            rmeConfig.shuffleMode=shuffleMode;
            rmeConfig.samPro=samPro;

            if  strcmpi(shuffleMode, 'noise')
                rmeConfig.randiHoldR=randiHoldR;
            else
                rmeConfig.randiHoldR=randiHold;
            end

            readmeFileName=strcat(outputFolderName,'readme.txt');
            done=writeReadmeHypo(readmeFileName, fixedTypes, rmeConfig);



            %% Node Analyze


            xcoordsTotal=G.Nodes.Coordinates(:, 1);
            ycoordsTotal=G.Nodes.Coordinates(:, 2);
            zcoordsTotal=G.Nodes.Coordinates(:, 3);

            %
            % figure
            % h=plot(G,'XData',xcoordsTotal,'YData',ycoordsTotal);



            intersectThr=0.25;

            pNodesCell=cell(numExtMotifs, 1);
            pHNodesCell=cell(numExtMotifs, 1);
            ppHNodesCell=cell(numExtMotifs, 1);
            ppHNodesUCell=cell(numExtMotifs, 1);


            ppHMNodesCell=cell(numExtMotifs, 1);



            ppHNodesCorr=zeros(numExtMotifs, 1);
            numelVec=zeros(numExtMotifs, 1);
            seedsToPWMC=cell(numExtMotifs, 1);
            % extract motif nodes

            for imt=1:numExtMotifs
                extMotifi=extMotif{imt};
                pNodesCell{imt}=(extMotifi.pNodes(:)).';
                pHNodesCell{imt}=(extMotifi.pHNodes(:)).';

                pNodesPPH=[pNodesCell{imt}, pHNodesCell{imt}];

                ppHNodesUCell{imt}=unique(pNodesPPH);
                numelVec(imt)=numel(ppHNodesUCell{imt});
                pNodesPPH=repelem(pNodesPPH, 1, ceil(exp(log(numExtMotifs-extMotifi.testPvalue))));

                ppHNodesCell{imt}=pNodesPPH;

                ppHMNodesCell{imt}=[extMotifi.pNodes;extMotifi.pHNodes];

                seedsToPWMC{imt}=extMotifi.seedsToPWM;


            end




            ppHNodesCorrMat=zeros(numExtMotifs);

            numNodesVec=zeros(numExtMotifs,1);

            % intersect between all motif pairs
            all_marks = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};

            figure

            for  refi=1:numExtMotifs
                ppHNodesCorr=zeros(numExtMotifs, 1);

                for iAll=1:numExtMotifs

                    ppHNodesCorr(iAll)=sum(ismember(ppHNodesUCell{iAll},ppHNodesUCell{refi}))/max(numelVec(iAll), numelVec(refi));
                end


                ppHNodesCorrMat(:, refi)=ppHNodesCorr;

                numNodesVec(refi)=numel(ppHNodesCell{iAll});


                ppHNodesCorr=sort(ppHNodesCorr, 'descend');
                plot(ppHNodesCorr,'Marker', all_marks{mod(refi-1,13)+1})
                %     text(2, ppHNodesCorr(2),num2str(refi))
                hold on

            end
            grid on

            xlim([1,numExtMotifs])

            legend(num2str((1:numExtMotifs).'), 'Location','eastoutside')
            xlabel('motif index')
            ylabel('node correlation')
            title('sorted node correlations')


            % figname=strcat('output\nodeIntersect', num2str(W), fignamExtnd);
            % saveas(gcf,figname)


            ppHNodesCorrMat=round(ppHNodesCorrMat, 2);

            figure

            %         heatmap(ppHNodesCorrMat);
            imagesc(ppHNodesCorrMat);
            title('motif node correlations heatmap')

            colormap(parula	)
            colorbar
            % strcat(outputFolderName,'evalueWG', num2str(W), fignamExtnd);
            figname=strcat(outputFolderName,'nodeIntersectHMap',fignamExtnd);
            saveas(gcf,figname)


                saveFileName=strcat(outputFolderName,'matData', filenamExtnd);
                save(saveFileName, '-regexp', '^(?!(extMotifC|hfig|lengthList|listLength|lengthListR|listLengthR|nodeListR|nHoldSeq|negSeq|pHoldSeq|posSeq|seqData|ppHNodesCell)$).')


            %% Highlight Intersections
            % green	[0 1 0]
            % magenta	[1 0 1]
            % yellow	[1 1 0]

            % red		[1 0 0]

            % black	[0 0 0]
            % cyan	[0 1 1]
            % white	[1 1 1]


            % if cgOptions.isRandom
            %     figure
            % else
            %     figure('WindowState' ,'maximize')
            % end



            for iMG=1:numExtMotifs

                figure('WindowState' ,'maximize')
                % figure

                ax1=subplot(4, 3, [1,2,4,5,7,8]); % top row of the 3x3 grid



                % figure('WindowState' ,'maximize')
                % figure


                % h3=plot(ax1,G,'XData',xcoordsTotal,'YData',ycoordsTotal, 'ZData',zcoordsTotal,'NodeColor',[0, 0.8, 0.8], 'EdgeColor',[0.8, 0.8, 0.8]);
                h3=plot(ax1,G,'XData',xcoordsTotal,'YData',ycoordsTotal, 'ZData',zcoordsTotal,'NodeColor',[0, 0.8, 0.8], 'EdgeColor','none');

                motifPaths=ppHMNodesCell{iMG};
                cMTypes=cPTypes(motifPaths);

                neronFlags=all(cMTypes>0, 2);
                if any(neronFlags)
                    cMTypes=cMTypes(neronFlags, :);
                end


                [cMTypesU, iu, ju]=unique(cMTypes, 'rows');
                ncolor=[1 0 1;0 1 0;1 0.64 0;1 0 0;0 0 0;0 1 1;rand(length(iu),3)];

                cCounts=accumarray(ju, 1);
                [cCounts, ist]=sort(cCounts, 'descend');
                cCountsTotal=sum(cCounts);
                cMTypesU=cMTypesU(ist, :);

                cli=ncolor(1, :);

                for iCTU=1:length(iu)
                    cli=ncolor(iCTU, :);
                    seedPaths=motifPaths(ju==ist(iCTU), :);
                    seedPaths=seedPaths(:);


                    % highlight(h3,seedPaths,'NodeColor',cli, 'EdgeColor',cli);
                    highlight(h3,seedPaths,'NodeColor',cli);

                end


                titleStr=sprintf("motif %d- total Counts %d", iMG, cCountsTotal);
                for iStr=1:min(length(cCounts), W+1)
                    countStr=sprintf(":%d ", cCounts(iStr));
                    titleStr=strcat(titleStr, countStr);
                end
                title(titleStr)




                hold on
                numLegend=min(6, length(iu));
                ax=zeros(numLegend, 1);

                for iLg=1:numLegend
                    ax(iLg)=plot(NaN,NaN,'.', 'MarkerFaceColor', ncolor(iLg, :), 'MarkerEdgeColor',  ncolor(iLg, :), 'markersize', 20); %plotting invisible points of desired colors
                end
                legendText=cTypeChars(cMTypesU(1:numLegend, :));
                legend(ax, legendText)




                hold off
                xlabel(ax1, 'x');
                ylabel(ax1, 'y');
                zlabel(ax1, 'z');


                ax2=subplot(4, 3, 3);
                ax3 = subplot(4, 3, 6);
                ax4 = subplot(4, 3, 9);
                ax5 = subplot(4, 3, 10);
                ax6 = subplot(4, 3, 11);
                ax7 = subplot(4, 3, 12);



                copyobj(h3,ax2);

                copyobj(h3,ax3);
                copyobj(h3,ax4);
                copyobj(h3,ax5);
                copyobj(h3,ax6);
                copyobj(h3,ax7);

                xlabel(ax2, 'x');
                ylabel(ax2, 'y');
                zlabel(ax2, 'z');

                xlabel(ax3, 'x');
                ylabel(ax3, 'y');
                zlabel(ax3, 'z');

                xlabel(ax4, 'x');
                ylabel(ax4, 'y');
                zlabel(ax4, 'z');

                xlabel(ax5, 'x');
                ylabel(ax5, 'y');
                zlabel(ax5, 'z');

                xlabel(ax6, 'x');
                ylabel(ax6, 'y');
                zlabel(ax6, 'z');


                xlabel(ax7, 'x');
                ylabel(ax7, 'y');
                zlabel(ax7, 'z');


                view(ax2,90, 0); title(ax2,'view along -x');
                % set(ax2, 'Xdir', 'reverse')

                view(ax3,270, 0); title(ax3,'view along +x');
                % set(ax3, 'Ydir', 'reverse')

                view(ax4,0, 0); title(ax4,'view along +y');

                view(ax5,0, 270); title(ax5,'view along +z');

                set(ax5, 'Ydir', 'reverse')




                view(ax6,0, 90); title(ax6,'view along -z');

                view(ax7,0, 180); title(ax7,'view along -y');
                set(ax7, 'Zdir', 'reverse')


                %
                % set(gca,'xtick',BregmaAllU*1e5/2,'xticklabel',BregmaAllU)
                % if length(BregmaAllU)<2
                %     set(gca,'ytick',(0:length(AnimalIDIdx)-1)*3000,'yticklabel',(1:length(AnimalIDIdx)))
                % else
                %     set(gca,'ytick',(AnimalIDIdx-1)*3000,'yticklabel',(1:length(AnimalIDIdx)))
                % end
                % xlabel('Bregma')
                % ylabel('Animal\_ID')

                figname=strcat(outputFolderName,'motifGroupOnGraph', num2str(iMG),num2str(W), fignamExtnd);
                % print(gcf,figname(1:end-4), '-djpeg', '-r600'); %<-Save as jpg with 600 DPI
                % print(gcf,figname(1:end-4), '-djpeg'); %<-Save as jpg

                saveas(gcf,figname)
            end




            check=1;

             %% Highlight Intersections All
            % green	[0 1 0]
            % magenta	[1 0 1]
            % yellow	[1 1 0]

            % red		[1 0 0]

            % black	[0 0 0]
            % cyan	[0 1 1]
            % white	[1 1 1]


            % if cgOptions.isRandom
            %     figure
            % else
            %     figure('WindowState' ,'maximize')
            % end
                figure('WindowState' ,'maximize')


                ax1=subplot(4, 3, [1,2,4,5,7,8]); % top row of the 3x3 grid
                ax2=subplot(4, 3, 3);
                ax3 = subplot(4, 3, 6);
                ax4 = subplot(4, 3, 9);
                ax5 = subplot(4, 3, 10);
                ax6 = subplot(4, 3, 11);
                ax7 = subplot(4, 3, 12);


                % figure('WindowState' ,'maximize')
                % figure


                % h3=plot(ax1,G,'XData',xcoordsTotal,'YData',ycoordsTotal, 'ZData',zcoordsTotal,'NodeColor',[0, 0.8, 0.8], 'EdgeColor',[0.8, 0.8, 0.8]);
                h3=plot(ax1,G,'XData',xcoordsTotal,'YData',ycoordsTotal, 'ZData',zcoordsTotal,'NodeColor',[0, 0.8, 0.8], 'EdgeColor','none');
                % ncolorAll=[0 1 1;1 0 1;0 1 0;1 0.64 0;1 0 0;0 0 0;rand(numExtMotifs,3)];
                ncolorAll=[0 0.4470 0.7410;0.9290 0.6940 0.1250;0.8500 0.3250 0.0980;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;0 0 0;0 1 1;1 0 1;0 1 0;1 0.64 0;1 0 0;rand(numExtMotifs,3)];

            for iMG=1:numExtMotifs

                % figure


                motifPaths=ppHMNodesCell{iMG};
                cMTypes=cPTypes(motifPaths);

                neronFlags=all(cMTypes>0, 2);
                if any(neronFlags)
                    cMTypes=cMTypes(neronFlags, :);
                end


                [cMTypesU, iu, ju]=unique(cMTypes, 'rows');

                ncolori=ncolorAll(iMG, :);


                ncolor=repmat(ncolori, length(iu), 1);

                cCounts=accumarray(ju, 1);
                % [cCounts, ist]=sort(cCounts, 'descend');
                cCountsTotal=sum(cCounts);

                cli=ncolor(1, :);

                for iCTU=1:length(iu)
                    cli=ncolor(iCTU, :);
                    seedPaths=motifPaths(ju==iCTU, :);
                    seedPaths=seedPaths(:);


                    % highlight(h3,seedPaths,'NodeColor',cli, 'EdgeColor',cli);
                    highlight(h3,seedPaths,'NodeColor',cli);

                end




            end


                % title("All Motifs Highlighted")

                hold off
                xlabel(ax1, 'x');
                ylabel(ax1, 'y');
                zlabel(ax1, 'z');




                copyobj(h3,ax2);

                copyobj(h3,ax3);
                copyobj(h3,ax4);
                copyobj(h3,ax5);
                copyobj(h3,ax6);
                copyobj(h3,ax7);

                xlabel(ax2, 'x');
                ylabel(ax2, 'y');
                zlabel(ax2, 'z');

                xlabel(ax3, 'x');
                ylabel(ax3, 'y');
                zlabel(ax3, 'z');

                xlabel(ax4, 'x');
                ylabel(ax4, 'y');
                zlabel(ax4, 'z');

                xlabel(ax5, 'x');
                ylabel(ax5, 'y');
                zlabel(ax5, 'z');

                xlabel(ax6, 'x');
                ylabel(ax6, 'y');
                zlabel(ax6, 'z');


                xlabel(ax7, 'x');
                ylabel(ax7, 'y');
                zlabel(ax7, 'z');


                view(ax2,90, 0); title(ax2,'view along -x');
                % set(ax2, 'Xdir', 'reverse')

                view(ax3,270, 0); title(ax3,'view along +x');
                % set(ax3, 'Ydir', 'reverse')

                view(ax4,0, 0); title(ax4,'view along +y');

                view(ax5,0, 270); title(ax5,'view along +z');

                set(ax5, 'Ydir', 'reverse')




                view(ax6,0, 90); title(ax6,'view along -z');

                view(ax7,0, 180); title(ax7,'view along -y');
                set(ax7, 'Zdir', 'reverse')


                %
                % set(gca,'xtick',BregmaAllU*1e5/2,'xticklabel',BregmaAllU)
                % if length(BregmaAllU)<2
                %     set(gca,'ytick',(0:length(AnimalIDIdx)-1)*3000,'yticklabel',(1:length(AnimalIDIdx)))
                % else
                %     set(gca,'ytick',(AnimalIDIdx-1)*3000,'yticklabel',(1:length(AnimalIDIdx)))
                % end
                % xlabel('Bregma')
                % ylabel('Animal\_ID')


             figname=strcat(outputFolderName,'allMotifGroupOnGraph', num2str(iMG),num2str(W), fignamExtnd);
                % print(gcf,figname(1:end-4), '-djpeg', '-r600'); %<-Save as jpg with 600 DPI
                % print(gcf,figname(1:end-4), '-djpeg'); %<-Save as jpg

               saveas(gcf,figname)


            check=1;

        end
    end
end


check=1;


% hmo.title='cell type involved min motifs with their percentage';
%% Beta Analysis
% clc
% close all

% cOTypes=cellTable.cellSubtypeOrig;
%
% xyCoords=[cellTable.Centroid_X, cellTable.Centroid_Y];
%
%
%
% for iMt=1:numExtMotifs
%
%     fprintf('evaluating motif number:  %d\n', iMt);
%
%     MotifPaths=ppHMNodesCell{iMt};
%
%
%
%     for iPost=1:W
%
%         MPTypesi=cPTypes(MotifPaths(:,iPost));
%
%         ANodeFlags=MPTypesi==cgOptions.ambigType;
%
%         if any(ANodeFlags)
%
%
%
%             [cntPost, MPTypesiU]=cntNodeCTypes(MotifPaths(:,iPost),cPTypes);
%
%
%
%             fig = figure('WindowState' ,'maximize');
%
%
%
%             subplot(1, 3, 1)
%
%             bar((1:length(cntPost)), cntPost)
%
%             xticks(1:length(cntPost))
%
%             set(gca,'XTickLabel' , num2cell(MPTypesiU));
%             ylim([0, sum(cntPost)])
%             yticks(unique([cntPost; sum(cntPost)]))
%             grid on
%             title(sprintf('motif position cell types (numNode=%d)',size(MotifPaths, 1)) )
%
%             cntPostPredict=cntPost(MPTypesiU~=cgOptions.ambigType);
%
%             cntPostPredict=round(cntPostPredict/sum(cntPostPredict)*100);
%
%             MPTypesiPrU=MPTypesiU(MPTypesiU~=cgOptions.ambigType);
%
%
%             [cntODif, AMPTypesOU]=cntNodeCTypes(MotifPaths(ANodeFlags,iPost),cOTypes);
%
%             subplot(1, 3, 2)
%
%             bar((1:length(cntODif)), cntODif)
%             xticks(1:length(AMPTypesOU))
%             ylim([0, sum(cntODif)])
%             yticks(unique([cntODif; sum(cntODif)]))
%             grid on
%
%
%             set(gca,'XTickLabel' , num2cell(AMPTypesOU));
%
%             title(sprintf('Ambig real types(numA=%d)', sum(ANodeFlags)))
%
%             subplot(1, 3, 3)
%
%             xyCoordsi=xyCoords(MotifPaths(ANodeFlags,iPost), :);
%
%             nearestNeighs=knnsearch(xyCoords, xyCoordsi,'K',2);
%
%             nearestNeighs=nearestNeighs(:, 2);
%
%             [cntNNPost, MPTypesNNiU]=cntNodeCTypes(nearestNeighs,cOTypes);
%
%
%             bar((1:length(cntNNPost)), cntNNPost)
%             xticks(1:length(cntNNPost))
%             cntNNPost=sort(cntNNPost);
%             set(gca,'XTickLabel' , num2cell(MPTypesNNiU));
%             ylim([0, sum(cntNNPost)])
%             yticks(unique([cntNNPost; sum(cntNNPost)]))
%             grid on
%
%             title(sprintf('Ambig nearest neighbor types(numA=%d)', sum(ANodeFlags)))
%
%             xlabelTex=sprintf('cell type, mt=%d, post=%d',iMt, iPost);
%
%             sam=axes(fig,'visible','off');
%             sam.XLabel.Visible='on';
%             sam.YLabel.Visible='on';
%             ylabel(sam,'frequency');
%             xlabel(sam,xlabelTex);
%
%             figname=strcat(outputFolderName,'ambigProbs', num2str(iMt),num2str(iPost), fignamExtnd);
%             % print(gcf,figname(1:end-4), '-djpeg', '-r600'); %<-Save as jpg with 600 DPI
%             % print(gcf,figname(1:end-4), '-djpeg'); %<-Save as jpg
%
%             saveas(gcf,figname)
%
%
%         end
%
%
%
%
%
%     end
%
%
%
% end



