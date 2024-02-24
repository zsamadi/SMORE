clc
close all
clearvars

if length(findall(0))>1
    delete(findall(0));
end


% v03: cretaed to match with hypo_v01, fig5_NRet
% v02: created for noise edition of control data generation
% the output results very much depend on the holdout and training dada, so
% be careful with that

cd(fileparts(which(mfilename)));
cd('..\')

addpath(genpath(pwd))

isTestMode=false;
isPlotGraph=true;

c=clock;
c=c(1:end-1);
cStr=string(c(1));
for ic=2:length(c)-1
    cStr=strcat(cStr, '_', string(c(ic)));
end

trNumShuffleo=1;




% delete(gcp('nocreate'))
% parpool;
rng(1650);
%% Hyperparameters
gTrainNum=50;
W=4;
% Hyperparameters for Creating Graph
RBPDiscard=false;           % Discard RBP Cells from downward analysis.

% Hyperparameters for Random Walk
numExtMotifs=20;
isEraseNodes=true;
fixedTypes=0;
% fixedTypes=[(0:96),(111:116)];

fixedTypes=fixedTypes(:);


%% Creat Graph

% cNodes=(1:1000).';
% ucNodes=(1:100).';
%
% unqMersFlags=getZNICIL(cNodes,  ones(W, 1),ucNodes);



folderName='D:\nucla\P1\matIO\data\';
plotGraph=true;

cgOptions.isJClusterID=false;
cgOptions.bregID=12;

cgOptions.ambigRatio=0;

cgOptions.ambigType=13;

cgOptions.isNoiseCTS=true;
cgOptions.fixedTypes=fixedTypes;

cgOptions.retinaAndSection=[3, 3];
cgOptions.ORetIDs=40;

cgOptions.plotGraph=false;
cgOptions.isUniformWeight=true;
cgOptions.isRandom=false;
cgOptions.nCTypesRand=12;
cgOptions.nNodesRand=cgOptions.nCTypesRand*840;
cgOptions.randCoords=false;
cgOptions.gMode="delaunay"; % "delaunay" or "knn", or "epsilon"


cgOptions.isNBreg=false;
cgOptions.isPBreg=true;
cgOptions.isFemale=true;
cgOptions.isMale=false;
cgOptions.isNaive=true;
cgOptions.isStimu=true;



outputFolderName0='D:\nucla\P1\matIO\output\bipolar\';

outputFolderName0=strcat(outputFolderName0,'W', num2str(W));

if ~isTestMode
    outputFolderName0=strcat(outputFolderName0,'_',cStr);
end





% outputFolderName=strcat(outputFolderName,'ID',num2str(cgOptions.ORetIDs));
%
% if cgOptions.isFemale && cgOptions.isMale
%     outputFolderName=strcat(outputFolderName, 'AseFM');
% elseif cgOptions.isFemale
%     outputFolderName=strcat(outputFolderName, 'AseF');
% else
%     outputFolderName=strcat(outputFolderName, 'AseM');
% end
%
%
% if cgOptions.isNBreg && cgOptions.isPBreg
%     outputFolderName=strcat(outputFolderName, 'BrgPN');
% elseif cgOptions.isNBreg
%     outputFolderName=strcat(outputFolderName, 'BrgN');
% else
%     outputFolderName=strcat(outputFolderName, 'BrgP');
% end
%
% if cgOptions.isNaive && cgOptions.isStimu
%     outputFolderName=strcat(outputFolderName, 'BhvNS');
% elseif cgOptions.isNaive
%     outputFolderName=strcat(outputFolderName, 'BhvN');
% else
%     outputFolderName=strcat(outputFolderName, 'BhvS');
% end
%
% outputFolderName=strcat(outputFolderName, cgOptions.gMode);
% outputFolderName=strcat(outputFolderName, shuffleMode);




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
    fignamExtnd='bipolar.jpeg';
    filenamExtnd='bipolar.mat';
    matFigExtend='bipolar.fig';
    textExtend='bipolar.txt';
    xGLLim=-6000;
    xGHLim=6000;
    yGLLim=-5500;
    yGHLim=3500;
end

% rmeConfig.fixedTypes=fixedTypes;
% rmeConfig.gMode=cgOptions.gMode;
% rmeConfig.shuffleMode=shuffleMode;

% done=writeReadmebipolar(outputFolderName,  rmeConfig);

cgOptions.rEpsilon=300;
cgOptions.xyNoiseStd=0;

iDumb=0;
gHLoptions.isPlotEdges=false;

kNeighsMax=15;

pvalueCell=cell(30,1);

cellTypesOne=["1a","1b","2","3a","3b","4","5a","5b","5c","5d","6","7","8","9","RBC" ];
cellTypesOne=cellTypesOne(:);


gHLoptions.is3D=false;
gHLoptions.ctAnnot=cellTypesOne;

kNeighs=5;

for shuffleMode="shuffle"%["shuffle", "kernelPath"] %can be 'shuffle', 'cluster', 'kernel', 'noisy'
    for fixedTypes=[0, 15]

         close all
         if length(findall(0))>1
             delete(findall(0));
          end

        cgOptions.numNeighs=kNeighs;
        [G,gStruct, ~]=creatGraph(folderName,cgOptions);

        trNumShuffle=trNumShuffleo;
        iDumb=iDumb+1;
        if iDumb<1
            continue
        end

        outputFolderName=strcat(outputFolderName0, '\', num2str(iDumb-1));


        outputFolderName=strcat(outputFolderName, '\');

        if ~exist(outputFolderName, 'dir')
            mkdir(outputFolderName)
        end


        cPTypes=(G.Nodes.label(:, 1)).';
        cSections=(G.Nodes.label(:, 2)).';
        cSectionsU=unique(cSections);
        xyCoordinates=G.Nodes.Coordinates;




        shOptions.rEpsilon=200;
        shOptions.fixedTypes=fixedTypes;
        shOptions.isSectShuffle=true;
        shOptions.numClusters=20;

        shConfig.cSections=cSections;
        shConfig.fixedTypes=fixedTypes;
        shConfig.xyCoordinates=xyCoordinates;


        if isPlotGraph
            gHLoptions.folderName=outputFolderName;
            gHLoptions.isShuffled=false;

            plotHLightGraph(G, cPTypes, gHLoptions);
        end



        %% Type Shuffling

        shConfig.fixedNodes=[];
        shConfig.shuffleMode=shuffleMode;

        if strcmpi(shuffleMode, 'noise')
            cgOptions.xyNoiseStd=3*edgeDistStd;
            [GR,gStructR]=creatGraph(folderName,cgOptions);

            cNTypes=(GR.Nodes.label(:, 1)).';
            cSectionsR=(GR.Nodes.label(:, 2)).';
            cSectionsUR=unique(cSectionsR);
            xyCoordinatesR=GR.Nodes.Coordinates;
            %
            %  notSampleVertex=(1:length(cPTypes));
            %
            % notSampleVertex= notSampleVertex(cPTypes==15);
            notSampleVertex=[];
            [nodeListR, lengthListR, nodeSectionsR]=enumerateUSRKPathsTest(gStructR, pathLength, [ones(1,pathLength-1), samPro], notSampleVertex);
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
            shConfig.nearNeighs=[];
            shConfig.numClusters=1;
            shConfig.numShuffle=trNumShuffleo;
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


        % for sectioni=cSectionsU
        %     sectionFlagi=(cSections==sectioni);
        %     cPTypesi=cPTypes(sectionFlagi);
        %     PWMT=sum(cPTypesi.'==cPTypesU);
        %     PWMT=PWMT/sum(PWMT);
        %
        %     figure
        %     bar(cPTypesU, PWMT)
        %     grid on
        %     xlabel('cell type')
        %     ylabel('frequency')
        %     text(cPTypesU-0.5, PWMT+0.001, cellTypesOne)
        %     title(sprintf('animal ID %d, Bregma ID',sectioni) )
        % end









        % figname='output\cellTypesFreqs.fig';
        %     % print(gcf,figname(1:end-4), '-djpeg', '-r600'); %<-Save as jpg with 600 DPI
        % print(gcf,figname(1:end-4), '-djpeg'); %<-Save as jpg
        %
        % saveas(gcf,figname)

        % cPTypes=cPTypes(randperm(length(cPTypes)));
        alphabet=cellTypesOne;
        numCells=length(alphabet);


        % highlighTypes(G, cPTypes, 'Primary')
        % highlighTypes(G, cNTypesCluster, 'Cluster Shuffled')
        % highlighTypes(G, cNTypesKernel, 'Kernel Shuffled')

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
        samPro0=0;

        pathLength=W;
        samPro=1;

        if samPro~=samPro0
            % notSampleVertex=(1:length(cPTypes));
            % notSampleVertex= notSampleVertex(cPTypes==15);
            notSampleVertex=[];
            enOptions.pdv=[ones(1,pathLength-1), samPro];
            enOptions.k=pathLength;
            enOptions.isChRadial=true;


            [pathList, lengthList, nodeSections, pathWeights]=enumerateUSRKPathsTest(gStruct,enOptions);

        end
        samPro0=samPro;
        %% Generate Shuffled data with k-mers presereved

        timerValue=tic;
        gOptions.hFrac=0;
        gOptions.mkvOrder=0;
        gOptions.rvp=true;
        % gOptions.isChar=false;
        % gOptions.cntFold=0; % used to generate larger number of control sequences, 2 fold, 3 fold, silly :)

        gOptions.isSectHold=true;

        gOptions.isHExclsv=true;

        % nodeSectionsN=nodeSections;


        gOptions.numNodes=length(cPTypes);

        nodeSectionsGen=nodeSections;
        nodeSectionsGenu=unique(nodeSectionsGen);

        randiHold=nodeSectionsGenu(randperm(length(nodeSectionsGenu)));
        randiHold=randiHold(1:floor(length(nodeSectionsGenu)/4));


        nodeSectionsGen(ismember(nodeSectionsGen,randiHold))=1;
        gOptions.numGRs=1;


        [posSeq, pHoldSeq, cPTypes0, posWeight, pHoldWeight]=generateSeqs(pathList,pathWeights,nodeSectionsGen,cPTypes, gOptions);
        % pHoldSeq=posSeq;
        shuffleModeGMaps=shuffleMode;


        if strcmpi(shuffleMode, 'noise')
            trNumShuffle=1;
            nodeSectionsGenR=nodeSectionsR;
            nodeSectionsGenRu=unique(nodeSectionsGenR);

            % randiHoldR=randperm(length(nodeSectionsGenRu));
            % randiHoldR=randiHoldR(1:floor(length(nodeSectionsGenRu)/4));
            randiHoldR=randiHold;
            nodeSectionsGenR(ismember(nodeSectionsGenR,randiHoldR))=1;

            [negSeq, nHoldSeq, cNTypes0]=generateSeqs(nodeListR,nodeSectionsGenR,cNTypes, gOptions);
            shuffleModeGMaps= "kernelPath";
        else
            gOptions.numGRs=trNumShuffle;
            [~, ~, cNTypes0, negWeight]=generateSeqs(pathList,pathWeights,nodeSectionsGen,cNTypes, gOptions);

            % if kNeighs==kNeighsMax
            %     [negSeq, ~, cNTypes0]=generateSeqs(nodeList,nodeSectionsGen,cNTypes, gOptions);
            % end

            % nHoldSeq=pHoldSeq;

        end

        delete(gcp('nocreate'))




        % gOptions.isControl=true;
        %
        % [negSeq, nHoldSeq]=generateSeqs(pWordList,nodeSectionsGen,gOptions);

        %% Identify Motifs



        % parpool;


        % alphabetNumeric=alphabet;



        % numDStream=size(posSeq, 1);


        % Write Fasta file to be analyzed by official STREME
        % numDStreamFasta=numDStream;
        % posFasta=char(posSeq+64);
        % negFasta=char(negSeq+64);

        % numDStreamFasta=5;
        %     rng('default');


        % alphabet=char(alphabetNumeric+64);
        % numHNodes=sum(gStruct.labels(:, 2)==1);



        % cNTypes=cPTypes(randperm(length(cPTypes)));
        indSeedMode=false;
        diffMotif=false;
        isUBack=false;

        [extMotif,textOut, commandText, background, extMotifEval]=mtStreme(cPTypes=cPTypes0, cNTypes=cNTypes0,cPHTypes=cPTypes,cNHTypes=cNTypes(1:length(cPTypes)),pSeq=posSeq, ...
            posWeight=posWeight, negWeight=negWeight, rvp=gOptions.rvp, mkvOrder=gOptions.mkvOrder, wMin=W, wMax=W,threshold=0.005, nmotifs=numExtMotifs,alphabet=alphabet, ...
            isEraseNodes=isEraseNodes,fixedTypes=fixedTypes,shConfig=shConfig,shuffleMode=shuffleModeGMaps,  nRefIter=20, isUBack=isUBack, trNumShuffle=trNumShuffle, diffMotif=diffMotif, indSeedMode=indSeedMode,...
             gTrainNum=gTrainNum, isNormalTest=false);


        elapsedTime=toc(timerValue);


        %%  resolve motif nodes
                seqData.pSeq=posSeq;
                seqData.nSeq=posSeq;
                seqData.pHSeq=posSeq;
                seqData.nHSeq=posSeq;
                seqData.back=background;
                seqData.cPTypes=cPTypes;
                seqData.cNTypes=cNTypes;
                seqData.cPHTypes=cPTypes;
                seqData.cNHTypes=cNTypes;
                seqData.rvp=gOptions.rvp;

                erzOptions.isErzNegative=false;

                erzOptions.fixedTypes=fixedTypes;

                erzOptions.shuffleMode=shuffleMode;
                erzOptions.nearNeighs=shConfig.nearNeighs;
                erzOptions.isErzFixNodes=false;
                erzOptions.isOLess=false;
                erzOptions.trNumShuffle=trNumShuffle;
                erzOptions.cPTypesInit=cPTypes;
                erzOptions.cNTypesInit=cNTypes;
                erzOptions.isErzHOut=false;
                erzOptions.isOKSingleNC=true;
                erzOptions.isNotNMer=true;
                erzOptions.isEraseNodes=true;
                erzOptions.rvp=gOptions.rvp;
                erzOptions.mkvOrder=gOptions.mkvOrder;


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
        filename=strcat(outputFolderName, 'zSTREME','_Hypo', num2str(W), '.txt');
        writeTextOutput(filename, extMotif,bkg,  textOut,commandText, elapsedTime);


        numExtMotifs=length(extMotif);

        ProfDateNum=datetime("now");
        numFixedTypes=length(fixedTypes);

        fprintf('shuffling method %s and numFixedTypes %d finished at : %s\n', shuffleMode, numFixedTypes, ProfDateNum);








        %% Seqlogo



        for iexMotif=1:numExtMotifs

            exMotifi=extMotif{iexMotif};
            exMotifi.background=background;

            [WL, hfig, cTypeChars] = seqlogoRet(exMotifi);
            [~, exMotifi.cSeed]=max(exMotifi.PWM);

            cSeedStr = strjoin(string(exMotifi.cSeed), '_');

            figname=strcat(outputFolderName, 'mLogo',num2str(iexMotif), '_',cSeedStr, fignamExtnd);
            saveas(hfig,figname)

        end

        %% PValues

        pvalueVec=zeros(numExtMotifs, 1);
        scoreThrVec=zeros(numExtMotifs, 1);

        for iexMotif=1:numExtMotifs

            exMotifi=extMotif{iexMotif};
            pvalueVec(iexMotif)=exMotifi.testPvalue;
            scoreThrVec(iexMotif)=exMotifi.scoreThr;

        end

        figure
        plot(pvalueVec, '-kv', 'LineWidth',1)
        grid on

        xlabel('motif index')
        ylabel('evalue')

        figname=strcat(outputFolderName,'evalueWG', num2str(W), fignamExtnd);
        saveas(gcf,figname)

        pvalueCell{iDumb}=pvalueVec;


        %
        % figure
        % plot(sort(scoreThrVec))
        % grid on

        rmeConfig=cgOptions;

        rmeConfig.cTypeChars=cTypeChars;

        rmeConfig.cellTypes=cellTypesOne;

        rmeConfig.W=W;

        rmeConfig.randiHold=randiHold;

        if  strcmpi(shuffleMode, 'noisy')
            rmeConfig.randiHoldR=randiHoldR;
        else
            rmeConfig.randiHoldR=randiHold;
        end

        readmeFileName=strcat(outputFolderName,'readme.txt');
        rmeConfig.isEraseFixed=false;
        done=writeReadmeNRet(readmeFileName, fixedTypes, rmeConfig);



        %% Node Analyze


        xcoordsTotal=G.Nodes.Coordinates(:, 1);
        ycoordsTotal=G.Nodes.Coordinates(:, 2);
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
        % all_marks = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};


        saveFileName=strcat(outputFolderName,'matData', filenamExtnd);
        save(saveFileName, '-regexp', '^(?!(extMotifC|hfig|lengthList|listLength|lengthListR|listLengthR|nodeListR|nodeSectionsR|nHoldSeq|negSeq|pHoldSeq|posSeq|seqData)$).')


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


            h3=plot(G,'XData',xcoordsTotal,'YData',ycoordsTotal);

            motifPaths=ppHMNodesCell{iMG};
            cMTypes=cPTypes(motifPaths);

            neronFlags=all(cMTypes>14, 2);
            if any(neronFlags)
                cMTypes=cMTypes(neronFlags, :);
            end


            [cMTypesU, iu, ju]=unique(cMTypes, 'rows');
            cCounts=accumarray(ju, 1);
            [cCounts, ist]=sort(cCounts, 'descend');
            cMTypesU=cMTypesU(ist, :);
            cCountsTotal=sum(cCounts);
            ncolor=[1 0 1;0 1 0;1 0.64 0;1 0 0;0 0 0;0 1 1;rand(length(iu),3)];


            cli=ncolor(1, :);

            for iCTU=1:length(iu)
                cli=ncolor(iCTU, :);
                seedPaths=motifPaths(ju==ist(iCTU), :);
                seedPaths=seedPaths(:);


                highlight(h3,seedPaths,'NodeColor',cli, 'EdgeColor',cli);
            end


            titleStr=sprintf("motif %d- total Counts %d", iMG, cCountsTotal);
            for iStr=1:min(length(cCounts), W+1)
                countStr=sprintf(":%d ", cCounts(iStr));
                titleStr=strcat(titleStr, countStr);
            end
            title(titleStr)

            grid minor
            sectionU=(1:3);
            AnimalIDIdx=(1:3);


            set(gca,'xtick',(sectionU-1)*15000,'xticklabel',sectionU)
            set(gca,'ytick',(AnimalIDIdx-1)*10000,'yticklabel',(1:length(AnimalIDIdx)))


            xlabel('section')
            ylabel('Animal\_ID')
            numLegend=min(6, length(iu));
            hold on
            ax=zeros(numLegend, 1);

            for iLg=1:numLegend
                ax(iLg)=plot(NaN,NaN,'.', 'MarkerFaceColor', ncolor(iLg, :), 'MarkerEdgeColor',  ncolor(iLg, :), 'markersize', 20); %plotting invisible points of desired colors
            end
            legendText=cTypeChars(cMTypesU(1:numLegend, :));
            legend(ax, legendText)



            figname=strcat(outputFolderName,'motifGroupOnGraph', num2str(iMG),num2str(W), fignamExtnd);
            % print(gcf,figname(1:end-4), '-djpeg', '-r600'); %<-Save as jpg with 600 DPI
            % print(gcf,figname(1:end-4), '-djpeg'); %<-Save as jpg

            saveas(gcf,figname)

        end
    end
end

check=1;


