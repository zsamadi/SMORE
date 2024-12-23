function smore(dataFolder , options)

arguments

    dataFolder.input='sampleDataEb.csv';
    dataFolder.gexFilename='sampleDataEbGEx.csv';
    dataFolder.gAnotFilename='';
    dataFolder.output='output/sampleDataEb';
    options.nMotifs=5;
    options.WMotif=4;
    options.ND=0;
    options.gMode='delaunay'; % "delaunay" or "knn", or "epsilon"
    options.nNeighs=5;
    options.rEps=0;
    options.samplingFreq=1;
    options.shuffleMode='shuffle';
    options.fixedTypes=0;
    options.neighDepth=4;
    options.nTrain=50;
    options.nScore=10;
    options.isEnrich=true;
    options.diffMotif=false;
    options.nEval=25;
    options.nRefine=4;
    options.nRefineIter=20;
    options.doGEA=true;
    options.gePvalMin=5e-2;
    options.iSPGEx=false;
    options.geRandTest=false;
    options.isGEByTissue=false;
    options.disObvs=false;
end

rng(1650);

isRunURPENSmore=true;

cd(fileparts(which(mfilename)));
cd('..\')

addpath(genpath(pwd))

if isRunURPENSmore
    outputFolderName=strcat(dataFolder.output, '/');
    
    outputFolderName=strcat(outputFolderName, string(datetime('today')), '/');
    mNodefilename=strcat(outputFolderName, 'ppHMNodesCell.mat');

else

    outputFolderName=dataFolder.output;

    

    mNodefilename=strcat(outputFolderName, 'ppHMNodesCell.mat');
    load(mNodefilename);
end


clc

LogicalStr = {'false', 'true'};

commandText=sprintf(['smore(input="%s",...\n output="%s",...\n nMotifs=%d,... \n WMotif=%d,' ...
    ['...\n ND=%d,...\n gMode="%s",...\n nNeighs=%d,...\n rEps=%d,...\n samplingFreq=%1.2f,' ...
    '...\n shuffleMode="%s",...\n fixedTypes=%d,...\n neighDepth=%d,...\n nScore=%d,...\n nTrain=%d,' ...
    '...\n diffMotif=%s,...\n isEnrich=%s,...\n nEval=%d,...\n nRefine=%d ,...\n '], ...
                   'nRefineIter=%d,...\n doGEA=%s,...\n gePvalMin=%1.3f,...\n geRandTest=%s) '], dataFolder.input,outputFolderName,options.nMotifs, ...
                   options.WMotif,options.ND,options.gMode,options.nNeighs,options.rEps, ...
                   options.samplingFreq, options.shuffleMode,options.fixedTypes,options.neighDepth, ...
                   options.nScore, options.nTrain,LogicalStr{options.diffMotif+1},LogicalStr{options.isEnrich+1}, ...
                   options.nEval,options.nRefine, options.nRefineIter,LogicalStr{options.doGEA+1},options.gePvalMin,LogicalStr{options.geRandTest+1});



statOut=sprintf('Smore started with shuffling method "%s" and FixedTypes "%s" \n', options.shuffleMode, num2str(options.fixedTypes));

fprintf(statOut);




isPlotGraph=true;



%% Creat Graph
W=options.WMotif;
isEraseNodes=true;

folderName=dataFolder.input;
cgOptions.isJClusterID=false;


cgOptions.isNoiseCTS=true;

cgOptions.retinaAndSection=[3, 3];

cgOptions.plotGraph=false;
cgOptions.isUniformWeight=true;
cgOptions.isRandom=false;
cgOptions.nCTypesRand=12;
cgOptions.nNodesRand=cgOptions.nCTypesRand*840;
cgOptions.randCoords=false;
cgOptions.gMode=options.gMode; 
cgOptions.nNeighs=options.nNeighs;
cgOptions.rEps=options.rEps;




fignamExtnd='.jpeg';



cgOptions.rEpsilon=0;
cgOptions.xyNoiseStd=0;


shuffleMode=options.shuffleMode;
fixedTypes=options.fixedTypes;

cgOptions.nNeighs=options.nNeighs;
cgOptions.ND=options.ND;
cgOptions.SIDSel=[];


[G,gStruct,cidNOut, ocellTypesOne, haveZ, tissueIDs, TOut]=creatAGraph(folderName,cgOptions);

if haveZ
    gHLoptions.is3D=true;
else
    gHLoptions.is3D=false;
end

trNumShuffle=1;



if ~exist(outputFolderName, 'dir')
    mkdir(outputFolderName)
end



cPTypes=(G.Nodes.label(:, 1)).';
cSections=(G.Nodes.label(:, 2)).';
xyCoordinates=G.Nodes.Coordinates;


shConfig.cSections=cSections;
shConfig.fixedTypes=fixedTypes;
shConfig.xyCoordinates=xyCoordinates;


cTypeChars = ['ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz', char((198:198+length(ocellTypesOne)-53))];


[~, ~, jCPU]=unique(cPTypes);
cCPU=accumarray(jCPU, 1);
[~, iSCU]=sort(cCPU, 'descend');

alphabet=cTypeChars(1:length(ocellTypesOne));
alphabet(iSCU)=alphabet;

cellTypesOne=string(ocellTypesOne(:));

%% plot graph


gHLoptions.isPlotEdges=false;
gHLoptions.ctAnnot=cellTypesOne;
gHLoptions.tissueIDs=tissueIDs;

gHLoptions.iSideView=false;
gHLoptions.alphabet=alphabet;






if isPlotGraph
    gHLoptions.folderName=outputFolderName;
    gHLoptions.isShuffled=false;


gHLoptions.ncolorL= {'414141';
'BD9F6D';
'F6EB14';
'33D4FF';
'ED047A';
'3C419A';
'86C0C8';
'8282BE';
'69BD45';
'33FFF1'};

gHLoptions.iShowPlot=true;



    plotHLightGraph(G, G.Nodes.label(:, 1), gHLoptions);
end


if isRunURPENSmore


%% Sample Paths from the Graph, Rand_ESU, or Random Walk

statOut=sprintf('Sampling input graph with %d nodes and %d edges, with sampling frequency %1.2f\n',length(gStruct.labels(:, 1)),length(G.Edges.EndNodes),options.samplingFreq);
sampleSpecs=sprintf('Sampled input graph with %d nodes and %d edges, with sampling frequency %1.2f\n',length(gStruct.labels(:, 1)),length(G.Edges.EndNodes),options.samplingFreq);

fprintf(statOut);

pathLength=W;
samPro=options.samplingFreq;


enOptions.pdv=[ones(1,pathLength-1), samPro];
enOptions.k=pathLength;
enOptions.isChRadial=true;
enOptions.disObvs=options.disObvs;
URPENTimerValue=tic;
enOptions.ctypes=cPTypes(:);
enOptions.seedsAll=[];

if W>2
    [pathList, ~, nodeSections, pathWeights]=URPENH(gStruct,enOptions);
else
    pathList=G.Edges.EndNodes;
    nodeSections=gStruct.labels(pathList(:, 1), 2);
    pathWeights=G.Edges.Weight;
    pathList=num2cell(pathList, 2);

end
URPENElapsedTime=toc(URPENTimerValue);

statURPENG=sprintf('Generated %d samples of length %d in %d seconds\n',length(pathList),options.WMotif, ceil(URPENElapsedTime));
fprintf(statURPENG);


%% Type Shuffling



if strcmpi(options.shuffleMode, 'kernel')
    statOut='Generate kernels for kernel shuffling methd \n';
    fprintf(statOut);
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


shConfig.partitionInd=[];
shConfig.fixedNodes=[];
shConfig.shuffleMode=shuffleMode;
shConfig.rEpsilon=cgOptions.rEpsilon;
shConfig.fixedTypes=fixedTypes;
shConfig.isSectShuffle=true;
shConfig.numClusters=20;
shConfig.cSections=cSections;
shConfig.numClusters=1;
shConfig.numShuffle=1;
cNTypes=getShuffleTypes(cPTypes, shConfig);
if isPlotGraph
    gHLoptions.folderName=outputFolderName;
    gHLoptions.isShuffled=true;
    gHLoptions.iShowPlot=false;

    plotHLightGraph(G, cNTypes(1:length(cPTypes)), gHLoptions);
end

cPTypesU=unique(cPTypes);
PWMT=sum(cPTypes(:)==(cPTypesU(:)).');
PWMT=PWMT/sum(PWMT);

figure('visible','off');
subplot(1, 2, 1)
bar(cPTypesU, PWMT)
grid on
xlabel('cell type')
ylabel('frequency')
text(cPTypesU-0.5, PWMT+0.001, alphabet(:))
title('Primary background')

PWMTNK=sum(cNTypes(:)==cPTypesU);
PWMTNK=PWMTNK/sum(PWMTNK);
subplot(1, 2, 2)
bar(cPTypesU, PWMTNK)
grid on
xlabel('cell type')
ylabel('frequency')
text(cPTypesU-0.5, PWMT+0.001, alphabet(:))
title('shuffled background')

figname=strcat(outputFolderName,'primShuffleBg', fignamExtnd);

saveas(gcf,figname)
close(gcf);


%% Generate Shuffled data

timerValue=tic;
gOptions.hFrac=0;
gOptions.mkvOrder=0;
gOptions.rvp=true;
gOptions.isSectHold=true;
gOptions.isHExclsv=true;
gOptions.numNodes=length(cPTypes);
gOptions.disObvs=options.disObvs;
nodeSectionsGen=nodeSections;
nodeSectionsGenu=unique(nodeSectionsGen);
randiHold=nodeSectionsGenu(randperm(length(nodeSectionsGenu)));
randiHold=randiHold(1:floor(length(nodeSectionsGenu)/4));
nodeSectionsGen(ismember(nodeSectionsGen,randiHold))=1;
gOptions.numGRs=1;
[posSeq, ~, cPTypes0, posWeight]=generateSeqs(pathList,pathWeights,nodeSectionsGen,cPTypes, gOptions);
shuffleModeGMaps=shuffleMode;
gOptions.numGRs=trNumShuffle;
[~, ~, cNTypes0, negWeight]=generateSeqs(pathList,pathWeights,nodeSectionsGen,cNTypes, gOptions);
delete(gcp('nocreate'))


%% Identify Motifs
rng(1650);

indSeedMode=false;
isUBack=false;

[extMotif,textOut, ~, background]=mtSmore(cPTypes=cPTypes0, cNTypes=cNTypes0,cPHTypes=cPTypes,cNHTypes=cNTypes(1:length(cPTypes)),pSeq=posSeq, ...
    posWeight=posWeight, negWeight=negWeight, rvp=gOptions.rvp, mkvOrder=gOptions.mkvOrder, wMin=W, wMax=W,threshold=0.005, nmotifs=options.nMotifs,alphabet=alphabet, ...
    isEraseNodes=isEraseNodes,fixedTypes=fixedTypes,shConfig=shConfig,shuffleMode=shuffleModeGMaps,  nRefIter=options.nRefineIter, isUBack=isUBack, trNumShuffle=trNumShuffle, diffMotif=options.diffMotif, indSeedMode=indSeedMode,...
    gTrainNum=options.nTrain, isEnrich=options.isEnrich, scIterMax=options.nScore, patience=3,NREF=options.nRefine, NEVAL=options.nEval);
elapsedTime=toc(timerValue);

%%  resolve motif nodes

if samPro<1
    enOptions.pdv=ones(1,pathLength);
    seedsAll=arrayfun(@(i) extMotif{i}.seedsToPWM, 1:length(extMotif), 'uniformOutput', false);
    seedsAll=vertcat(seedsAll{:});

    enOptions.pdv=ones(1,pathLength);
    enOptions.seedsAll=seedsAll;
    enOptions.ctypes=cPTypes;

    if W>2
        [pathList, ~, nodeSections, pathWeights]=URPENH(gStruct,enOptions);
        
    else
        pathList=G.Edges.EndNodes;
        nodeSections=gStruct.labels(pathList(:, 1), 2);
        pathWeights=G.Edges.Weight;
        pathList=num2cell(pathList, 2);

    end
    nodeSectionsGen=nodeSections;
    nodeSectionsGen(ismember(nodeSectionsGen,randiHold))=1;
    posSeq=generateSeqs(pathList,pathWeights,nodeSectionsGen,cPTypes, gOptions);
end

seqData.pSeq=posSeq;
seqData.nSeq=[];
seqData.pHSeq=[];
seqData.nHSeq=[];
seqData.back=background;
seqData.cPTypes=cPTypes;
seqData.cNTypes=[];
seqData.cPHTypes=[];
seqData.cNHTypes=[];
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
filename=strcat(outputFolderName, 'outputResults.txt');
writeTextOutput(filename, extMotif,bkg,  textOut,commandText, elapsedTime, alphabet);


numExtMotifs=length(extMotif);

ProfDateNum=datetime("now");
numFixedTypes=sum(fixedTypes>0);

statOut=sprintf('Shuffling method "%s" and numFixedTypes %d finished at : %s\n', shuffleMode, numFixedTypes, ProfDateNum);
fprintf(statOut)

statOut=sprintf('logos and motif highlights are saved in "%s" \n', string(outputFolderName));
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
    close(hfig);

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
ylabel('pvalue')

figname=strcat(outputFolderName,'pvalueW', num2str(W), fignamExtnd);
saveas(gcf,figname)
close(gcf);



%

rmeConfig=cgOptions;

rmeConfig.cTypeChars=cTypeChars;

rmeConfig.cellTypes=(1:length(cellTypesOne));

rmeConfig.cellTypesOne=cellTypesOne;

rmeConfig.W=W;

rmeConfig.randiHold=randiHold;

if  strcmpi(shuffleMode, 'noisy')
    rmeConfig.randiHoldR=randiHoldR;
else
    rmeConfig.randiHoldR=randiHold;
end

readmeFileName=strcat(outputFolderName,'readme.txt');
rmeConfig.isEraseFixed=false;
rmeConfig.sampleSpecs=sampleSpecs;
rmeConfig.statURPENG=statURPENG;
writeReadme(readmeFileName, fixedTypes, rmeConfig);

%% Node Analyze



ppHMNodesCell=cell(numExtMotifs, 1);
seedsToPWMC=cell(numExtMotifs, 1);

for imt=1:numExtMotifs
    extMotifi=extMotif{imt};


    ppHMNodesCell{imt}=[extMotifi.pNodes;extMotifi.pHNodes];

    seedsToPWMC{imt}=extMotifi.seedsToPWM;

end

save(mNodefilename,'ppHMNodesCell', 'extMotif', 'options');

%% Highlight Motif Nodes

hlConfig.is3D=false;
hlConfig.W=W;

hlConfig.cTypeChars=cTypeChars;

hlConfig.fignamExtnd=fignamExtnd;


highLMNodes(G, ppHMNodesCell,outputFolderName, hlConfig)

if haveZ
    hlConfig.is3D=true;

    highLMNodes(G, ppHMNodesCell,outputFolderName, hlConfig)
end



%% Highlight all motifs 
hlConfig.is3D=false;
hlConfig.numHLights=min(10, numExtMotifs);
hlConfig.ndesPlotIDx=[];

highALMNodes(G, ppHMNodesCell,outputFolderName, hlConfig)

if haveZ
    hlConfig.is3D=true;

    highALMNodes(G, ppHMNodesCell,outputFolderName, hlConfig)
end

end


%% read cellTable
if options.doGEA
if options.iSPGEx
    geNames=importdata(dataFolder.gAnotFilename);

    mCells=length(cidNOut);
    nGenes=length(geNames);


    gexData = h5read(dataFolder.gexFilename, '/X/data');
    indCols = h5read(dataFolder.gexFilename, '/X/indices');
    indRPtr = h5read(dataFolder.gexFilename, '/X/indptr');

    indRows=repelem((1:mCells).', indRPtr(2:end)-indRPtr(1:end-1));


    gExpression=sparse(indRows, indCols-min(indCols)+1, gexData,mCells, nGenes);

    gExpression=gExpression(cidNOut, :);




else



    gExpression=readtable(dataFolder.gexFilename);
    geNames=gExpression.Properties.VariableNames;


    if strcmpi(geNames{1}, 'CID')
        geNames=geNames(2:end);
        gExpression=gExpression(:, 2:end);
    end

    gExpression=gExpression(cidNOut, :);
    gExpression=table2array(gExpression);






end
end


%% Gene Expression analysis


if abs(options.gePvalMin)>1
    pvalAllHMapMin=exp(-abs(options.gePvalMin));
else
    pvalAllHMapMin=options.gePvalMin;
end




if options.doGEA 


    statOut='Starting gene expression analysis \n';
    fprintf(statOut)
    geaSpecs.isPlotPDF=false;
    geaSpecs.isAddNoise=true;
    geaSpecs.pvalAllHMapMin=pvalAllHMapMin;
    geaSpecs.outputFolderName=strcat(outputFolderName, 'genEx/');
    geaSpecs.isTwoSided=true; 
    geaSpecs.isRndTest=options.geRandTest;
    geaSpecs.W=W;
    geaSpecs.geneAnnotes=geNames;
    geaSpecs.cellType=cPTypes(:);
    geaSpecs.isGEByTissue=options.isGEByTissue;

    geaSpecs.SID=G.Nodes.label(:, 2);
    geaSpecs.alphabet=alphabet;
  
    if ~exist(geaSpecs.outputFolderName, 'dir')
         mkdir(geaSpecs.outputFolderName)        
    end
    if isRunURPENSmore
        gexStart=tic;
    
        gExResults=genExSPAnalysis(gExpression,ppHMNodesCell, geaSpecs);
        gextiming=toc(gexStart);
        save(mNodefilename,'ppHMNodesCell', 'extMotif', 'options', 'gExResults');
    end
    %%
    geaSpecs.pvalHMapMin=100;
    geaSpecs.dMedHMapMin=1;
   
    
    gExHPlot(gExResults, geaSpecs);

end



