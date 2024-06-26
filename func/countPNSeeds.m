function [seedsWmax, unqMersNC, weightsM, seqData]=countPNSeeds(seqData, options)

options.trNumShuffle=1;
shConfig=options.shConfig;
shConfig.numShuffle=1;

options.isPos=true;
options.isPN=false;
options.isOnlyPScore=false;
[seedsWmax, seqData]=countHPSeeds(seqData, options);
options.isPos=false;
unqMersNC=cell(options.gTrainNum, 1);

numPSeeds=size(seedsWmax.seeds, 1);
weightsM=zeros(numPSeeds, options.gTrainNum);

addPos=options.addPos;

numSeeds=zeros(options.gTrainNum+addPos, 1);
NSeedsC=cell(options.gTrainNum+addPos, 1);
NWeightsC=cell(options.gTrainNum+addPos, 1);
NCountsC=cell(options.gTrainNum+addPos, 1);

options.back=seqData.back;

options.isNotNMer=false;


NSeedsC{1}=seedsWmax.seeds;

NWeightsC{1}=seedsWmax.weights(:, 1);
numSeeds(1)=size(seedsWmax.seeds, 1);
NCountsC{1}=seedsWmax.counts(:, 1);


for iter=1:options.gTrainNum

    fprintf('%3d%%', round(iter/options.gTrainNum*100));
    seqData.cNTypes=options.cNTypesMat(:, iter);
    options.pnWeights=seqData.negWeight;
    unqMersN=countHPSeeds(seqData, options);
    NSeedsC{iter+addPos}=unqMersN.seeds;
    NWeightsC{iter+addPos}=unqMersN.weights(:, 1);
    NCountsC{iter+addPos}=unqMersN.counts(:, 1);
    numSeeds(iter+addPos)=size(unqMersN.seeds, 1);

    fprintf('\b\b\b\b')


end

allSeeds=vertcat(NSeedsC{:});
allNWeights=vertcat(NWeightsC{:});
allCounts=vertcat(NCountsC{:});

[allSeedsNU, ~, jU]=unique(allSeeds, 'rows');

allCountsNMean=accumarray(jU, allCounts)/(options.gTrainNum+addPos);

allWeightsNMean=accumarray(jU, allNWeights)/(options.gTrainNum+addPos);

allWeightsNVar=accumarray(jU, allNWeights.^2)/(options.gTrainNum+addPos)-allWeightsNMean.^2;

allWeightsNVar=max(allWeightsNVar, 0);

allSeeds=[seedsWmax.seeds;allSeedsNU];

[allSeedsU, ~, ~]=unique(allSeeds, 'rows');

[allNFlag, iAllNSeed]=ismember(allSeedsU,allSeedsNU, 'rows');
iAllNSeed=iAllNSeed(allNFlag);

nWeights=zeros(size(allSeedsU, 1), 1);
nWeights(allNFlag)=allWeightsNMean(iAllNSeed)*(options.gTrainNum+addPos);

nCounts=zeros(size(allSeedsU, 1), 1);
nCounts(allNFlag)=allCountsNMean(iAllNSeed)*(options.gTrainNum+addPos);

nWSTDs=zeros(size(allSeedsU, 1), 1);
nWSTDs(allNFlag)=sqrt(allWeightsNVar(iAllNSeed));

[allPFlag, iAllPSeed]=ismember(allSeedsU,seedsWmax.seeds, 'rows');
iAllPSeed=iAllPSeed(allPFlag);



pWeights=zeros(size(allSeedsU, 1), 1);
pWeights(allPFlag)=seedsWmax.weights(iAllPSeed, 1);

pCounts=zeros(size(allSeedsU, 1), 1);
pCounts(allPFlag)=seedsWmax.counts(iAllPSeed, 1);

seedsWmax.seeds=allSeedsU;

seedsWmax.weights=[pWeights, nWeights];
seedsWmax.wSTDs=nWSTDs;
seedsWmax.counts=[pCounts, nCounts];

maxWidth=zeros(size(allSeedsU, 1), 1);
maxWidth(allPFlag)=seedsWmax.maxWidth;
seedsWmax.maxWidth=maxWidth;

D15=options.numCells.^(options.wMax-1:-1:0);

seedsIndex=(sum(allSeedsU.*D15,2));
seedsWmax.seedsIndex=seedsIndex;

seedsWmax={seedsWmax};
