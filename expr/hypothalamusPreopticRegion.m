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
trNumShuffleo=1;


W=4;
numExtMotifs=20;
gTrainNum=1;

cgOptions.isJClusterID=false;
cgOptions.gMode="delaunay"; % "delaunay" or "knn", or "epsilon"

iCodes=[{'0100'},{'00000101'},{'001011'},{'001010'},{'01110111'},{'00100'},{'0011000'},{'01111'},...
    {'0110010'},{'0000101'},{'000011001'},{'000011000'},{'0000011'},{'0000111'},{'00001101'},{'001101'},...
    {'0110011'},{'011000'},{'0011001'},{'01110110'},{'0011111'},{'00111100'},{'000000'},{'0000100'},{'0101'},...
    {'0111001'},{'00010'},{'0111000'},{'00111101'},{'0001100'},{'00000100'},{'0001101'},{'001110'},{'0111010'},{'000111'},{'011110'},{'01101'}];


iTypes=[{'I-1'}, {'I-2'},{'I-3'},{'I-4'},{'I-5'},{'I-6'},{'I-7'},{'I-8'},{'I-9'},...
    {'I-10'},{'I-11'}, {'I-12'},{'I-13'},{'I-14'},{'I-15'},{'I-16'},{'I-17'},{'I-18'},{'I-19'},...
    {'I-20'},{'I-21'}, {'I-22'},{'I-23'},{'I-24'},{'I-25'},{'I-26'},{'I-27'},{'I-29'},{'I-30'},{'I-31'},...
    {'I-32'}, {'I-33'},{'I-34'},{'I-35'},{'I-36'},{'I-37'},{'H-1'}];

eTypes=[{'E-1'}, {'E-2'},{'E-3'},{'E-4'},{'E-5'},{'E-6'},{'E-7'},{'E-8'},{'E-9'},...
    {'E-10'},{'E-11'}, {'E-12'},{'E-13'},{'E-14'},{'E-15'},{'E-16'},{'E-17'},{'E-18'},{'E-19'},...
    {'E-20'},{'E-21'}, {'E-22'},{'E-23'},{'E-24'},{'E-25'},{'E-26'},{'E-27'},{'E-28'},{'E-30'}];


eCodes=[{'1110101'},{'111001'},{'1110100'},{'110100'},{'10110000'},{'110011'},{'110101'},{'110111'},...
    {'10100'},{'101100011'},{'11000'},{'110110011'},{'1001'},{'10101'},{'110110010'},{'1111010'},{'1101101'},...
    {'111100'},{'1111011'},{'10111'},{'101001'},{'1011001'},{'1000'},{'11011000'},{'101100010'},{'111000'},{'111011'},{'11111'},{'101101'}];

ieTypes=[iTypes, eTypes];

ieCodes=[iCodes, eCodes];

shConfig.fixedNodes=[];
bregID=1;
iDumb=0;

% testMat=[0, 0, 0;0, 0, 1;0, 1, 0;0, 1, 1;1, 0, 0;1, 0, 1]>0;
% testMat=[0, 0, 0;0, 0, 1;0, 1, 0;0, 1, 1]>0;
% testMat=[0, 0, 0;0, 0, 1;0, 1, 0;0, 1, 1]>0;

testMat=[0, 0, 0;0, 0, 1]>0;
% testMat=[0, 0, 0]>0;

samPro=1;

outputFolderName0='D:\nucla\P1\matIO\output\hypo\';

outputFolderName0=strcat(outputFolderName0,'W', num2str(W));


if ~isTestMode
    c=clock;
    c=c(1:end-1);
    cStr=string(c(1));
    for ic=2:length(c)-1
        cStr=strcat(cStr, '_', string(c(ic)));
    end


    outputFolderName0=strcat(outputFolderName0,'_',cStr);
end


outputFolderName0=strcat(outputFolderName0, '\');
kNeighsMax=5;

pvalueCell=cell(30,1);
for testId=1:size(testMat, 1)
    for bregID=1
        for shuffleMode=["kernelPath", "shuffle"] % ["shuffle", "kernel","noise", "kernelPath"]
            

            for kNeighs=kNeighsMax:-1:5
                cgOptions.numNeighs=kNeighs;

                iDumb=iDumb+1;
                if iDumb<1
                    continue;
                end
                trNumShuffle=trNumShuffleo;
                close all
                if length(findall(0))>1
                    delete(findall(0));
                end

                ambigRatio=0;


                testi=testMat(testId, :);

                isOLess=false;

                isEnrich=true;

                % isEnrich=~isOLess && isEnrich;
                isErzFixNodes=0;

                numFixedTypes=14*double(~testi(3));

                rng(1750);
                %% Hyperparameters

                isEraseNodes=true;
                fixedTypes=(0:numFixedTypes);
                fixedTypes=fixedTypes(:);


                %% Creat Graph


                folderName='D:\nucla\P1\matIO\data\';
                plotGraph=true;
                cgOptions.bregID=bregID;

                cgOptions.ambigRatio=ambigRatio;

                cgOptions.ambigType=84;
                cgOptions.ambgSWTs=[26,57,31];
                cgOptions.swRatio=0;

                cgOptions.isNoiseCTS=false;
                cgOptions.fixedTypes=fixedTypes;

                cgOptions.retinaAndSection=[3, 3];
                cgOptions.hypoIDs=40;

                cgOptions.plotGraph=true;
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


                outputFolderName=strcat(outputFolderName0, num2str(iDumb-1), '\');


                if ~exist(outputFolderName, 'dir')
                    mkdir(outputFolderName)
                end

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
                    fignamExtnd='Hypo.jpeg';
                    filenamExtnd='Hypo.mat';
                    matFigExtend='Hypo.fig';
                    textExtend='Hypo.txt';
                    xGLLim=-6000;
                    xGHLim=6000;
                    yGLLim=-5500;
                    yGHLim=3500;
                end


                [T, cellTypesOut, ambigType]=readCellTable(folderName, cgOptions);
                cgOptions.ambigType=ambigType;
                cellTypesOne=vertcat(cellTypesOut(1:14, 2), cellTypesOut(15:end, 1));
                [Lia,Locb] = ismember(ieTypes, cellTypesOne);

                cellCodes=cell(length(cellTypesOne), 1);
                cellCodes(Locb)=ieCodes;


                cgOptions.xyNoiseStd=0;
                cgOptions.cellTypesOne=cellTypesOne;
                cgOptions.rEpsilon=300;

                [G,gStruct, ~, BregmaAllU, AnimalIDIdx, cellTable,edgeDistStd]=creatHPGraph(T,outputFolderName,cgOptions);


                if isNoiseTest
                    randPermi=randperm(length(G.Nodes.label(:, 1)));
                    G.Nodes.label(:, 1)=G.Nodes.label(randPermi, 1);
                end
                cPTypes=(G.Nodes.label(:, 1)).';

                cSections=(G.Nodes.label(:, 2)).';
                cSectionsU=unique(cSections);



                %% Sample Paths from the Graph, Rand_ESU, or Random Walk
                 pathLength=W;

                enuOptions.k=pathLength;
                enuOptions.pdv=[ones(1,pathLength-1), samPro];
                enuOptions.isChRadial=true;

                [pathList, lengthList, nodeSections, pathWeights]=enumerateUSRKPathsTest(gStruct, enuOptions);


%%

neighNodesC=cell(length(gStruct), 1);

pathListAll=vertcat(pathList{:});

neighDepth=2;


for v=1:length(gStruct.neighs)

    if ~any(cPTypes(v)==fixedTypes)
        neighNodesT=v;

        for iDepth=1:neighDepth
            neighNodesT=gStruct.neighs(neighNodesT);
            neighNodesT=vertcat(neighNodesT{:});
        end
        neighNodesT=unique(neighNodesT);


        neighNodesFlag=ismember(pathListAll(:), neighNodesT);
        neighNodesFlag=reshape(neighNodesFlag, size(pathListAll));
        neighNodesFlag=any(neighNodesFlag, 2);



    
        neighNodes=pathListAll(neighNodesFlag, :);
        neighNodes=unique(neighNodes(:));
        neighNodesTypes=cPTypes(neighNodes);
        neighNodes(ismember(neighNodesTypes, fixedTypes))=[];
    
        neighNodesC{v}=neighNodes;
    end





end

pathListAll=[];

shConfig.nearNeighs=neighNodesC;
neighNodesC=[];

                %% Type Shuffling



                shConfig.partitionInd=[];


                shConfig.shuffleMode=shuffleMode;

                    shConfig.rEpsilon=cgOptions.rEpsilon;
                    shConfig.fixedTypes=fixedTypes;
                    shConfig.isSectShuffle=true;
                    shConfig.numClusters=20;
                    shConfig.cSections=cSections;
                    shConfig.fixedTypes=fixedTypes;
                    shConfig.numClusters=1;
                    shConfig.numShuffle=trNumShuffle;

                    cNTypes=getShuffleTypes(cPTypes, shConfig);






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
                [posSeq, pHoldSeq, cPTypes0, posWeight, pHoldWeight]=generateSeqs(pathList,pathWeights,nodeSectionsGen,cPTypes, gOptions);

                shuffleModeGMaps=shuffleMode;



                    gOptions.numGRs=trNumShuffle;
                    [~, ~, cNTypes0, negWeight]=generateSeqs(pathList,pathWeights,nodeSectionsGen,cNTypes, gOptions);



                %% Identify Motifs



                indSeedMode=false;
                diffMotif=true;
                isUBack=false;

                [extMotif,textOut, commandText, background, extMotifEval]=mtStreme(cPTypes=cPTypes0, cNTypes=cNTypes0,cPHTypes=cPTypes,cNHTypes=cNTypes(1:length(cPTypes)),pSeq=posSeq, ...
                    posWeight=posWeight, negWeight=negWeight, rvp=gOptions.rvp, mkvOrder=gOptions.mkvOrder, wMin=W, wMax=W,threshold=0.005, nmotifs=numExtMotifs,alphabet=alphabet, ...
                    isEraseNodes=isEraseNodes,fixedTypes=fixedTypes,shConfig=shConfig,shuffleMode=shuffleModeGMaps,  nRefIter=20, isUBack=isUBack, trNumShuffle=trNumShuffle, diffMotif=diffMotif, indSeedMode=indSeedMode,...
                    scIterMax=2, gTrainNum=gTrainNum, isNormalTest=~cgOptions.isUniformWeight);




                elapsedTime=toc(timerValue);

    

                %%  resolve motif nodes

                % posSeq=vertcat(posSeq{:});
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

                    [WL, hfig, cTypeChars]=seqlogoHypo(exMotifi);

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

                if gOptions.hFrac>0
                    rmeConfig.randiHold=randiHold;
                else
                    rmeConfig.randiHold=0;
                end

                rmeConfig.shuffleMode=shuffleMode;

                if  strcmpi(shuffleMode, 'noise')
                    rmeConfig.randiHoldR=randiHoldR;
                else
                    rmeConfig.randiHoldR=randiHold;
                end
                rmeConfig.samPro=samPro;

                readmeFileName=strcat(outputFolderName,'readme.txt');
                done=writeReadmeHypo(readmeFileName, fixedTypes, rmeConfig);

                pvalueCell{iDumb}=pvalueVec;



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

                    % pNodesPPH=[pNodesCell{imt}, pHNodesCell{imt}];

                    % ppHNodesUCell{imt}=unique(pNodesPPH);
                    % numelVec(imt)=numel(ppHNodesUCell{imt});
                    % pNodesPPH=repelem(pNodesPPH, 1, ceil(exp(log(numExtMotifs-extMotifi.testPvalue))));

                    % ppHNodesCell{imt}=pNodesPPH;

                    ppHMNodesCell{imt}=extMotifi.pNodes;

                    seedsToPWMC{imt}=extMotifi.seedsToPWM;


                end




                ppHNodesCorrMat=zeros(numExtMotifs);

                numNodesVec=zeros(numExtMotifs,1);

                % intersect between all motif pairs
                all_marks = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};





                saveFileName=strcat(outputFolderName,'matData', filenamExtnd);
                save(saveFileName, '-regexp', '^(?!(extMotifC|hfig|lengthList|listLength|lengthListR|listLengthR|nodeListR|nHoldSeq|negSeq|pHoldSeq|posSeq|seqData)$).')


                %
                % xcoordsTotal=G.Nodes.Coordinates(:, 1);
                % ycoordsTotal=G.Nodes.Coordinates(:, 2);
                % %
                % % figure
                % % h=plot(G,'XData',xcoordsTotal,'YData',ycoordsTotal);
                %
                %
                %
                % intersectThr=0.25;
                %
                % pNodesCell=cell(numExtMotifs, 1);
                % pHNodesCell=cell(numExtMotifs, 1);
                % ppHNodesCell=cell(numExtMotifs, 1);
                % ppHNodesUCell=cell(numExtMotifs, 1);
                %
                %
                % ppHMNodesCell=cell(numExtMotifs, 1);
                %
                %
                %
                % ppHNodesCorr=zeros(numExtMotifs, 1);
                % numelVec=zeros(numExtMotifs, 1);
                % seedsToPWMC=cell(numExtMotifs, 1);
                % % extract motif nodes
                %
                % for imt=1:numExtMotifs
                %     extMotifi=extMotif{imt};
                %     pNodesCell{imt}=(extMotifi.pNodes(:)).';
                %     pHNodesCell{imt}=(extMotifi.pHNodes(:)).';
                %
                %     pNodesPPH=[pNodesCell{imt}, pHNodesCell{imt}];
                %
                %     ppHNodesUCell{imt}=unique(pNodesPPH);
                %     numelVec(imt)=numel(ppHNodesUCell{imt});
                %     pNodesPPH=repelem(pNodesPPH, 1, ceil(exp(log(numExtMotifs-extMotifi.testPvalue))));
                %
                %     ppHNodesCell{imt}=pNodesPPH;
                %
                %     ppHMNodesCell{imt}=[extMotifi.pNodes;extMotifi.pHNodes];
                %
                %     seedsToPWMC{imt}=extMotifi.seedsToPWM;
                %
                %
                % end
                %
                %
                %
                %
                % ppHNodesCorrMat=zeros(numExtMotifs);
                %
                % numNodesVec=zeros(numExtMotifs,1);
                %
                % % intersect between all motif pairs
                % all_marks = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};
                %
                % figure
                %
                % for  refi=1:numExtMotifs
                %     ppHNodesCorr=zeros(numExtMotifs, 1);
                %
                %     for iAll=1:numExtMotifs
                %
                %         ppHNodesCorr(iAll)=sum(ismember(ppHNodesUCell{iAll},ppHNodesUCell{refi}))/max(numelVec(iAll), numelVec(refi));
                %     end
                %
                %
                %     ppHNodesCorrMat(:, refi)=ppHNodesCorr;
                %
                %     numNodesVec(refi)=numel(ppHNodesCell{iAll});
                %
                %
                %     ppHNodesCorr=sort(ppHNodesCorr, 'descend');
                %     plot(ppHNodesCorr,'Marker', all_marks{mod(refi-1,13)+1})
                %     %     text(2, ppHNodesCorr(2),num2str(refi))
                %     hold on
                %
                % end
                % grid on
                %
                % xlim([1,numExtMotifs])
                %
                % legend(num2str((1:numExtMotifs).'), 'Location','eastoutside')
                % xlabel('motif index')
                % ylabel('node correlation')
                % title('sorted node correlations')
                %
                %
                % % figname=strcat('output\nodeIntersect', num2str(W), fignamExtnd);
                % % saveas(gcf,figname)
                %
                %
                % ppHNodesCorrMat=round(sqrt(ppHNodesCorrMat), 2);
                %
                % figure
                %
                %  heatmap(ppHNodesCorrMat,'XLabel' , 'motif index','YLabel' , 'motif index','FontSize',12,'Colormap',parula);
                % % imagesc(ppHNodesCorrMat);
                % title('motif node correlations heatmap')
                %
                % colormap(parula)
                % colorbar
                % % strcat(outputFolderName,'evalueWG', num2str(W), fignamExtnd);
                % figname=strcat(outputFolderName,'nodeIntersectHMap',fignamExtnd);
                % saveas(gcf,figname)
                %
                % GCorr=graph(double(ppHNodesCorrMat>0));
                % motifBins = conncomp(GCorr);
                % [~, srtIDX]=sort(motifBins);
                %
                % ppHNodesCorrMat=ppHNodesCorrMat(srtIDX, :);
                % ppHNodesCorrMat=ppHNodesCorrMat(:, srtIDX);
                %
                % figure
                %
                %  heatmap(ppHNodesCorrMat,'XLabel' , 'motif index','YLabel' , 'motif index','FontSize',12,'Colormap',parula);
                % Ax = gca;
                % Ax.XDisplayLabels = srtIDX;
                % Ax.YDisplayLabels = srtIDX;
                %
                % % imagesc(ppHNodesCorrMat);
                % title('motif node correlations heatmap')
                %
                % colormap(parula)
                % colorbar
                % % strcat(outputFolderName,'evalueWG', num2str(W), fignamExtnd);
                % figname=strcat(outputFolderName,'nodeIntersectHMapOrdered',fignamExtnd);
                % saveas(gcf,figname)





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

                    neronFlags=all(cMTypes>0, 2);
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

                    set(gca,'xtick',BregmaAllU*1e5/2,'xticklabel',BregmaAllU)
                    if length(BregmaAllU)<2
                        set(gca,'ytick',(0:length(AnimalIDIdx)-1)*3000,'yticklabel',(1:length(AnimalIDIdx)))
                    else
                        set(gca,'ytick',(AnimalIDIdx-1)*3000,'yticklabel',(1:length(AnimalIDIdx)))
                    end
                    xlabel('Bregma')
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


                check=1;

                %% Dendrogram Analysis
                % clc
                % close all

                % denDist=cell(numExtMotifs*W, 1);
                % denDSpec=cell(numExtMotifs*W, 1);
                % denDistM=zeros(numExtMotifs*W, 6);
                % denDSpecC1=cell(numExtMotifs*W, 1);
                % iDend=1;
                % simLarC=cell(numExtMotifs, 1);
                % celLarC=cell(numExtMotifs, 1);
                % simLarCC=cell(numExtMotifs, 1);

                % for iMt=1:numExtMotifs
                %
                %     fprintf('evaluating motif number:  %d\n', iMt);
                %
                %     exMotifi=extMotif{iMt};
                %
                %     seedsToPWMi=exMotifi.seedsToPWM;
                %     seedsToPWMiU=unique(seedsToPWMi);
                %     seedsToPWMiU=seedsToPWMiU(seedsToPWMiU>14);
                %
                %     [f, ipos]=ismember(seedsToPWMi(:), seedsToPWMiU);
                %
                %     iposv=repmat((1:W), size(seedsToPWMi, 1));
                %     iposv=iposv(:);



                %
                % if iMt==3
                %     check=1;
                % end
                % simLari=zeros(length(seedsToPWMiU));
                % simLariC=cell(length(seedsToPWMiU));
                % pCellCodes=cellCodes(seedsToPWMiU);


                %
                %     for ii1=1:length(seedsToPWMiU)
                %
                %         iposv1=iposv(ipos==ii1);
                %         iposv1=unique(iposv1);
                %
                %         simLariC{ii1, 1}=[seedsToPWMiU(ii1), iposv1.'];
                %
                %
                %         for ii2=ii1+1:length(seedsToPWMiU)
                %
                %
                %             cmpLen=min(length(pCellCodes{ii1}), length(pCellCodes{ii2}));
                %             cmpLMax=max(length(pCellCodes{ii1}), length(pCellCodes{ii2}));
                %
                %             iCmp=1;
                %             while(pCellCodes{ii1}(iCmp)==pCellCodes{ii2}(iCmp))
                %                 iCmp=iCmp+1;
                %             end
                %
                %             simLari(ii1, ii2)=(iCmp-2)/(cmpLen-1);
                %             simLariC{ii1, ii2}=[iCmp-2, cmpLen-1, cmpLMax-1];
                %
                %
                %         end
                %     end
                %
                %     simLarC{iMt}=simLari;
                %     celLarC{iMt}=seedsToPWMiU;
                %     simLarCC{iMt}=simLariC;
                %
                % end
                %
                % denDist=denDist(1:iDend-1, :);
                % denDSpec=denDSpec(1:iDend-1, :);
                %
                % denDistM=denDistM(1:iDend-1, :);
                % denDSpecC1=denDSpecC1(1:iDend-1, :);



                ICellU=["I-23";"I-32";"I-2";"I-13";"I-24";"I-10";"I-12";"I-11";"I-15";"I-14";"I-27";"I-31";"I-33";"I-36";"I-6";"I-4";"I-3";"I-7";"I-19";"I-16";"I-34";...
                    "I-22";"I-30";"I-21";"I-1";"I-25";"I-18";"I-9";"I-17";"H-1";"I-29";"I-26";"I-35";"I-20";"I-5";"I-37";"I-8"];

                ECellU=["E-23";"E-13";"E-9";"E-14";"E-5";"E-25";"E-10";"E-22";"E-30";"E-20";"E-11";"E-21";"E-6";"E-4";"E-7";"E-24";"E-15";"E-12";"E-17";"E-8";"E-26";...
                    "E-2";"E-3";"E-1";"E-27";"E-18";"E-16";"E-19";"E-28"];

                [~, idxDng]=ismember([ICellU;ECellU], cellTypesOne);

                allIndxs=(1:length(cellTypesOne));
                valIndxs=idxDng;

                numTypes=length([ICellU;ECellU]);
                mt=zeros(numExtMotifs, numTypes);
                iMtp=1;
                mtp=zeros(W*numExtMotifs, numTypes);
                labelp=zeros(W*numExtMotifs, 1);

                for iMt=1:numExtMotifs

                    fprintf('evaluating motif number:  %d\n', iMt);

                    MotifPaths=ppHMNodesCell{iMt};
                    MPTypesi=cPTypes(MotifPaths(:));
                    % MPTypesi=MPTypesi(MPTypesi>14);
                    [MPTypesiU, iU, jU]=unique(MPTypesi);
                    tCnts=accumarray(jU,1);

                    [flagi, idxi]=ismember(MPTypesiU,valIndxs);
                    idxi=idxi(flagi);
                    tCnts=tCnts(flagi).';
                    mt(iMt, idxi.')=round(100*tCnts/sum(tCnts));

                    for iPst=1:W
                        MotifPathsi=MotifPaths(:,iPst);
                        MPTypesi=cPTypes(MotifPathsi);
                        [MPTypespiU, iU, jU]=unique(MPTypesi);
                        tCnts=accumarray(jU,1);
                        [flagi, idxi]=ismember(MPTypespiU,valIndxs);
                        idxi=idxi(flagi);
                        tCnts=tCnts(flagi).';
                        if length(tCnts)>1
                            mtp(iMtp, idxi.')=round(100*tCnts/sum(tCnts));
                            labelp(iMtp)=iMt*10+iPst;

                            iMtp=iMtp+1;
                        end
                    end

                end

                mtp=mtp(1:iMtp-1, :);
                labelp=labelp(1:iMtp-1);



                dgLevels=[1, 1.5, 2, 2.5, 3, 3.5, 4, 4.25, 4.5].';







                % gEZSIMeans=(1:37).';
                % gEZSIMeans(2)=gEZSIMeans(2)+0.1;
                %
                % gEZSIMeans(4)=gEZSIMeans(4)-0.01;
                % gEZSEMeans=(1:29).';

                ZIMn=[2,3;5,6;7,8;12,13;16, 17;18,19;22, 23;25, 26;28,29;31,32;34,35;36,37;4,38;...
                    9,40;14,41;15,42;20,43;24,44;27,46;33,48;50,1;51,10;52,11;55,21;56,30;57,47;63,...
                    49;64,62;65,45;54,61;67,53;39,59;69,58;70,60;71,68;72,66];



                ZIMn(ZIMn>37)=ZIMn(ZIMn>37)+29;

                ZIMnD=[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.2,1.2,1.2,1.2,1.2,1.2,1.3,1.4,1.5,1.3,1.4,1.3,1.4,1.5,1.6,1.7].';




                ZEMn=[1, 2;3, 4;6, 7;32, 5;33, 8;34, 9;35, 10;36, 31;37,30;12,13;11,39;14,15;...
                    17,18;16,42;19,43;20,44;45,41;46,40;21,22;23,24;25,49;50,48;27,28;52,26;29,53;54,51;55,47;56,38];
                ZEMn=ZEMn+37;

                ZEMn(ZEMn>66)=ZEMn(ZEMn>66)+73-37;

                ZEMnD=[1, 1, 1, 1.1,1.2,1.3,1.4, 1.5, 1.6, 1, 1.1, 1, 1, 1.1, 1.2,1.3,1.4,1.5,1,1.1,1.2,1.4,1.1,1.2,1.3,1.4,1.5,1.7].';



                ZMIE=[ZIMn;ZEMn;102,130];



                % ZI = linkage(gEZSIMeans,'complete');
                % D = pdist(gEZSIMeans);
                %
                % IleafOrder = optimalleaforder(ZI,D);
                %
                %
                % ZE=linkage(gEZSEMeans,'complete');
                %
                % ZIE=linkage([gEZSIMeans+100;gEZSEMeans],'complete');
                %
                % ZIE(end)=max(3, max(ZIE(1:end-1, end)));
                % ZIE(:, end)=log(1+ZIE(:, end));

                ZIE=zeros(size(ZMIE, 1),3);
                ZIE(:, 1:2)=ZMIE;

                dists=[ZIMnD;ZEMnD;1.7];

                [distsU, ~, jUD]=unique(dists);

                distsU=dgLevels(1:length(distsU));

                dists=distsU(jUD);

                ZIE(:, end)=dists;

                % figure()
                % dendrogram(ZI,0,  'Labels',ICellU)



                % dendrogram(ZI, 0)

                %
                % figure()
                % dendrogram(ZE, 0, 'Labels',ECellU)
                % dendrogram(ZE, 0)
                % xlim([0.5 inf])
                check=1;
                fig = figure('WindowState' ,'maximize');
                t=tiledlayout(8,1,'TileSpacing','tight', 'Padding','tight');
                nexttile([1,1])
                dendrogram(ZIE,0,  'Labels',[ICellU;ECellU], 'Reorder', (1:numTypes))
                xlim([0.5 numTypes+0.5])
                ylim([0 max(dists)])
                title('cell type involved in motifs with their percentage')

                nexttile([7, 1])
                hmo=heatmap(mt,'XLabel' , 'cell types','FontSize',12,'Colormap',cool);


                Ax = gca;
                Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
                Ax.XDisplayLabels=["   ";nan(length(Ax.XDisplayData)-1,1)];

                linkaxes(t.Children,'x')

                hmo.YLabel = 'extracted motifs';
                figname=strcat(outputFolderName,'cellTHeat', fignamExtnd);
                saveas(gcf,figname)




                fig = figure('WindowState' ,'maximize');
                t=tiledlayout(8,1,'TileSpacing','tight', 'Padding','tight');
                nexttile([1,1])
                dendrogram(ZIE,0,  'Labels',[ICellU;ECellU], 'Reorder', (1:numTypes))
                xlim([0.5 numTypes+0.5])
                ylim([0 max(dists)])
                title('cell type involved in motifs with their percentage')

                nexttile([7, 1])
                [mt, sidx]=sortrows(mt, 'descend');
                hmo=heatmap(mt,'XLabel' , 'cell types','FontSize',12,'Colormap',cool);


                Ax = gca;
                Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
                Ax.XDisplayLabels=["   ";nan(length(Ax.XDisplayData)-1,1)];
                Ax.YDisplayLabels=sidx;

                linkaxes(t.Children,'x')

                hmo.YLabel = 'extracted motifs';
                figname=strcat(outputFolderName,'cellTHeatOrdered', fignamExtnd);
                saveas(gcf,figname)



                fig = figure('WindowState' ,'maximize');
                t=tiledlayout(8,1,'TileSpacing','tight', 'Padding','tight');
                nexttile([1,1])
                dendrogram(ZIE,0,  'Labels',[ICellU;ECellU], 'Reorder', (1:numTypes))
                xlim([0.5 numTypes+0.5])
                ylim([0 max(dists)])
                title('cell type involved in motifs with their percentage')

                nexttile([7, 1])
                hmo=heatmap(mtp,'XLabel' , 'cell types','FontSize',12,'Colormap',cool);


                Ax = gca;
                Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
                Ax.XDisplayLabels=["   ";nan(length(Ax.XDisplayData)-1,1)];
                Ax.YDisplayLabels = labelp;

                linkaxes(t.Children,'x')

                hmo.YLabel = 'extracted motifs';
                figname=strcat(outputFolderName,'cellTHeatPost', fignamExtnd);
                saveas(gcf,figname)


                % hmo.title='cell type involved min motifs with their percentage';
                %% Beta Analysis
                % clc
                % close all

                cOTypes=cellTable.cellSubtypeOrig;

                xyCoords=[cellTable.Centroid_X, cellTable.Centroid_Y];



                for iMt=1:numExtMotifs

                    fprintf('evaluating motif number:  %d\n', iMt);

                    MotifPaths=ppHMNodesCell{iMt};



                    for iPost=1:W

                        MotifPathsi=MotifPaths(:,iPost);
                        MotifPathsi=unique(MotifPathsi);

                        MPTypesi=cPTypes(MotifPathsi);

                        ANodeFlags=MPTypesi==cgOptions.ambigType;

                        if any(ANodeFlags)



                            [cntPost, MPTypesiU]=cntNodeCTypes(MotifPathsi,cPTypes);



                            fig = figure('WindowState' ,'maximize');



                            subplot(1, 3, 1)

                            bar((1:length(cntPost)), cntPost)

                            xticks(1:length(cntPost))

                            set(gca,'XTickLabel' , num2cell(MPTypesiU));
                            ylim([0, sum(cntPost)])
                            yticks(unique([cntPost; sum(cntPost)]))
                            grid on
                            title(sprintf('motif position cell types (numNode=%d)',length(MotifPathsi)) )

                            cntPostPredict=cntPost(MPTypesiU~=cgOptions.ambigType);

                            cntPostPredict=round(cntPostPredict/sum(cntPostPredict)*100);

                            MPTypesiPrU=MPTypesiU(MPTypesiU~=cgOptions.ambigType);


                            [cntODif, AMPTypesOU]=cntNodeCTypes(MotifPathsi(ANodeFlags),cOTypes);

                            subplot(1, 3, 2)

                            bar((1:length(cntODif)), cntODif)
                            xticks(1:length(AMPTypesOU))
                            ylim([0, sum(cntODif)])
                            yticks(unique([cntODif; sum(cntODif)]))
                            grid on


                            set(gca,'XTickLabel' , num2cell(AMPTypesOU));

                            title(sprintf('Ambig real types(numA=%d)', sum(ANodeFlags)))

                            subplot(1, 3, 3)

                            xyCoordsi=xyCoords(MotifPathsi(ANodeFlags), :);

                            nearestNeighs=knnsearch(xyCoords, xyCoordsi,'K',2);

                            nearestNeighs=nearestNeighs(:, 2);

                            [cntNNPost, MPTypesNNiU]=cntNodeCTypes(nearestNeighs,cOTypes);


                            bar((1:length(cntNNPost)), cntNNPost)
                            xticks(1:length(cntNNPost))
                            cntNNPost=sort(cntNNPost);
                            set(gca,'XTickLabel' , num2cell(MPTypesNNiU));
                            ylim([0, sum(cntNNPost)])
                            yticks(unique([cntNNPost; sum(cntNNPost)]))
                            grid on

                            title(sprintf('Ambig nearest neighbor types(numA=%d)', sum(ANodeFlags)))

                            xlabelTex=sprintf('cell type, mt=%d, post=%d',iMt, iPost);

                            sam=axes(fig,'visible','off');
                            sam.XLabel.Visible='on';
                            sam.YLabel.Visible='on';
                            ylabel(sam,'frequency');
                            xlabel(sam,xlabelTex);

                            figname=strcat(outputFolderName,'ambigProbs', num2str(iMt),num2str(iPost), fignamExtnd);
                            % print(gcf,figname(1:end-4), '-djpeg', '-r600'); %<-Save as jpg with 600 DPI
                            % print(gcf,figname(1:end-4), '-djpeg'); %<-Save as jpg

                            saveas(gcf,figname)


                        end





                    end



                end

            end
        end
    end
end
