function smore(dataFolder , options)

arguments
    dataFolder.input='sampleDataEb.csv';
    dataFolder.output='output';
    options.alphabet=[];
    options.nEval=25;
    options.nRefine=4;
    options.WMotif=4;
    options.nRefineIter=20;
    options.threshold=0.005;
    options.patience=3;
    options.nMotifs=5;
    options.isErzHOut=true;
    options.nScore=1;
    options.nTrain=50;
    options.fixedTypes=0;
    options.diffMotif=false;
    options.shuffleMode='shuffle';
    options.isEnrich=true;
    options.samplingFreq=1;
    options.gMode='delaunay';
    options.numNeighs=5;
    options.rEps=0;
    options.neighDepth=4;
    options.ND=0;

end
clc

LogicalStr = {'false', 'true'};

commandText=sprintf(['smore(input="%s",...\n output="%s",... \n fixedTypes=%d,... \n  nEval=%d,... \n nRefine=%d ,... \n WMotif=%d,... \n ', ...
                   'nRefineIter=%d,...\n nmotifs=%d,... \n nScore=%d,... \n  nTrain=%d,...\n ', ...
                   'diffMotif=%s,...\n shuffleMode="%s",...\n isEnrich=%s,... \n samplingFreq=%1.2f,...\n gMode="%s",... \n numNeighs=%d,... ', ...
                   '\n rEps=%d,... \n neighDepth=%d)'], dataFolder.input,dataFolder.output,options.fixedTypes,options.nEval,options.nRefine, ...
                   options.WMotif,options.nRefineIter,options.nMotifs,options.nScore, ...
                   options.nTrain,LogicalStr{options.diffMotif+1}, options.shuffleMode,LogicalStr{options.isEnrich+1},options.samplingFreq, ...
                   options.gMode,options.numNeighs,options.rEps,options.neighDepth);



statOut=sprintf('smore started with shuffling method "%s" and numFixedTypes "%d" \n', options.shuffleMode, options.fixedTypes);

fprintf(statOut);

if strcmpi(options.shuffleMode, 'kernel')
    options.shuffleMode='kernelPath';
end

cd(fileparts(which(mfilename)));
cd('..\')

addpath(genpath(pwd))

isPlotGraph=true;


trNumShuffleo=1;


%% Hyperparameters
W=options.WMotif;
isEraseNodes=true;
fixedTypes=0;

fixedTypes=fixedTypes(:);


%% Creat Graph


folderName=dataFolder.input;
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
cgOptions.rEps=options.rEps;

outputFolderName=strcat(dataFolder.output, '\');



fignamExtnd='.jpeg';



cgOptions.rEpsilon=0;
cgOptions.xyNoiseStd=0;


shuffleMode=options.shuffleMode;
fixedTypes=options.fixedTypes;

% close all
% if length(findall(0))>1
%     delete(findall(0));
% end

cgOptions.numNeighs=options.numNeighs;
cgOptions.ND=options.ND;
[G,gStruct,~, nodeTypes]=creatAGraph(folderName,cgOptions);

trNumShuffle=trNumShuffleo;


% outputFolderName=strcat(outputFolderName0, '\');

if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName)
end


cPTypes=(G.Nodes.label(:, 1)).';
cSections=(G.Nodes.label(:, 2)).';
xyCoordinates=G.Nodes.Coordinates;


shConfig.cSections=cSections;
shConfig.fixedTypes=fixedTypes;
shConfig.xyCoordinates=xyCoordinates;

gHLoptions.isPlotEdges=false;

if isempty(options.alphabet)
    cTypeChars = ['ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz', char((198:198+length(nodeTypes)-53))];
    cellTypesOne=cTypeChars(1:length(nodeTypes));
else
    cellTypesOne=options.alphabet;
end
alphabet=cellTypesOne;
cellTypesOne=string(cellTypesOne(:));


gHLoptions.is3D=false;
gHLoptions.ctAnnot=cellTypesOne;


if isPlotGraph
    gHLoptions.folderName=outputFolderName;
    gHLoptions.isShuffled=false;

    plotHLightGraph(G, cPTypes, gHLoptions);
end






%% Sample Paths from the Graph, Rand_ESU, or Random Walk

pathLength=W;
samPro=options.samplingFreq;


enOptions.pdv=[ones(1,pathLength-1), samPro];
enOptions.k=pathLength;
enOptions.isChRadial=true;


[pathList, ~, nodeSections, pathWeights]=URPEN(gStruct,enOptions);

%%

if strcmpi(options.shuffleMode, 'kernel')
    shConfig.nearNeighs=cell(length(gStruct), 1);    
    pathListAll=vertcat(pathList{:});    
    neighDepth=options.neighDepth;    
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
            shConfig.nearNeighs{v}=neighNodes;
        end
    
    end
    pathListAll=[];
else
    shConfig.nearNeighs=[];
end

%% Type Shuffling
shConfig.partitionInd=[];
shConfig.fixedNodes=[];
shConfig.shuffleMode=shuffleMode;
shConfig.rEpsilon=cgOptions.rEpsilon;
shConfig.fixedTypes=fixedTypes;
shConfig.isSectShuffle=true;
shConfig.numClusters=20;
shConfig.cSections=cSections;
shConfig.fixedTypes=fixedTypes;
shConfig.numClusters=1;
shConfig.numShuffle=trNumShuffleo;
cNTypes=getShuffleTypes(cPTypes, shConfig);
if isPlotGraph
    gHLoptions.folderName=outputFolderName;
    gHLoptions.isShuffled=true;
    plotHLightGraph(G, cNTypes(1:length(cPTypes)), gHLoptions);
end

cPTypesU=(1:length(cellTypesOne));
PWMT=sum(cPTypes(:)==(cPTypesU(:)).');
PWMT=PWMT/sum(PWMT);

figure('visible','off');
subplot(1, 2, 1)
bar(cPTypesU, PWMT)
grid on
xlabel('cell type')
ylabel('frequency')
text(cPTypesU-0.5, PWMT+0.001, cellTypesOne)
title('Primary background')
% figname=strcat(outputFolderName,'PrimaryBg', fignamExtnd);

PWMTNK=sum(cNTypes(:)==cPTypesU);
PWMTNK=PWMTNK/sum(PWMTNK);
subplot(1, 2, 2)
bar(cPTypesU, PWMTNK)
grid on
xlabel('cell type')
ylabel('frequency')
text(cPTypesU-0.5, PWMT+0.001, cellTypesOne)
title('kernel shuffled background')

figname=strcat(outputFolderName,'kernelShfBg', fignamExtnd);

saveas(gcf,figname)


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
randiHold=nodeSectionsGenu(randperm(length(nodeSectionsGenu)));
randiHold=randiHold(1:floor(length(nodeSectionsGenu)/4));
nodeSectionsGen(ismember(nodeSectionsGen,randiHold))=1;
gOptions.numGRs=1;
[posSeq, ~, cPTypes0, posWeight]=generateSeqs(pathList,pathWeights,nodeSectionsGen,cPTypes, gOptions);shuffleModeGMaps=shuffleMode;
gOptions.numGRs=trNumShuffle;
[~, ~, cNTypes0, negWeight]=generateSeqs(pathList,pathWeights,nodeSectionsGen,cNTypes, gOptions);
delete(gcp('nocreate'))


%% Identify Motifs
rng(1650);

indSeedMode=false;
isUBack=false;

[extMotif,textOut, ~, background]=mtSmore(cPTypes=cPTypes0, cNTypes=cNTypes0,cPHTypes=cPTypes,cNHTypes=cNTypes(1:length(cPTypes)),pSeq=posSeq, ...
    posWeight=posWeight, negWeight=negWeight, rvp=gOptions.rvp, mkvOrder=gOptions.mkvOrder, wMin=W, wMax=W,threshold=options.threshold, nmotifs=options.nMotifs,alphabet=alphabet, ...
    isEraseNodes=isEraseNodes,fixedTypes=fixedTypes,shConfig=shConfig,shuffleMode=shuffleModeGMaps,  nRefIter=options.nRefineIter, isUBack=isUBack, trNumShuffle=trNumShuffle, diffMotif=options.diffMotif, indSeedMode=indSeedMode,...
    gTrainNum=options.nTrain, isEnrich=options.isEnrich, scIterMax=options.nScore, patience=options.patience,NREF=options.nRefine, NEVAL=options.nEval);
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
    extMotif{iEx}=outMotif;

end


bkg.PWM0=background{1};
bkg.mkvOrder=gOptions.mkvOrder;
filename=strcat(outputFolderName, 'out.txt');
writeTextOutput(filename, extMotif,bkg,  textOut,commandText, elapsedTime, alphabet);


numExtMotifs=length(extMotif);

ProfDateNum=datetime("now");
numFixedTypes=length(fixedTypes);

statOut=sprintf('shuffling method %s and numFixedTypes %d finished at : %s\n', shuffleMode, numFixedTypes, ProfDateNum);
fprintf(statOut)


%% Seqlogo

cTypeChars=alphabet;

for iexMotif=1:numExtMotifs

    exMotifi=extMotif{iexMotif};
    exMotifi.background=background;

    [~, hfig] = seqlogoGen(exMotifi, alphabet);
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

figure('visible','off');
plot(pvalueVec, '-kv', 'LineWidth',1)
grid on

xlabel('motif index')
ylabel('evalue')

figname=strcat(outputFolderName,'evalueWG', num2str(W), fignamExtnd);
saveas(gcf,figname)



%

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
writeReadme(readmeFileName, fixedTypes, rmeConfig);

%% Node Analyze


xcoordsTotal=G.Nodes.Coordinates(:, 1);
ycoordsTotal=G.Nodes.Coordinates(:, 2);
pNodesCell=cell(numExtMotifs, 1);
pHNodesCell=cell(numExtMotifs, 1);
ppHNodesCell=cell(numExtMotifs, 1);
ppHNodesUCell=cell(numExtMotifs, 1);
ppHMNodesCell=cell(numExtMotifs, 1);
numelVec=zeros(numExtMotifs, 1);
seedsToPWMC=cell(numExtMotifs, 1);
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




% saveFileName=strcat(outputFolderName,'matData', filenamExtnd);
% save(saveFileName, '-regexp', '^(?!(extMotifC|hfig|lengthList|listLength|lengthListR|listLengthR|nodeListR|nodeSectionsR|nHoldSeq|negSeq|pHoldSeq|posSeq|seqData)$).')


%% Highlight Intersections


for iMG=1:numExtMotifs


    figure('visible','off');

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

    saveas(gcf,figname)

end




