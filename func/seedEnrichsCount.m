function [zoopsCnt,totalCnt] =seedEnrichsCount(seedsWMax,psSeedIndex)

% from initial NEVAL motifs, chooses NREF motifs by counting all matches
% and approximate matches in primary and control data to the motifs and computing their pvalue
% input arguments:
%   - XnMerZoops: Cell containing zoops WMers of primary and control
%   sequences
%   - appxSeeds: Nm Seeds of length W that shoud be counted inside the zoops sites
%   - initial: true or false,  if true, only total counts are calculated,
%   if false, accumulated counts are calculated
%  Output arguments:
%   - seedCount: zoops count of the appxSeeds in the input WMers

numSeeds=length(psSeedIndex);


[~, indx]=ismember(seedsWMax.pSeedsIndex, psSeedIndex);
seedsSites=cell(numSeeds, 1);
seedsSites(indx)=seedsWMax.pSiteCell;
seedsWeights=cell(numSeeds, 1);

seedsWeights(indx)=seedsWMax.pWeightCell;


pWeights=cell(numSeeds, 1);
pWeights(indx)=seedsWMax.pWeightCell;
pWeights=vertcat(pWeights{:});


seedsCounts=zeros(length(seedsSites), 1);

for iSt=1:length(seedsWeights)
    seedsCounts(iSt)=sum(seedsWeights{iSt});
end


W=size(seedsWMax.seeds, 2);

allSeedsIndex=seedsWMax.seedsIndex;


allSeedIndexg=gpuArray(allSeedsIndex);
psSeedIndexg=gpuArray(psSeedIndex);
[~,seedStFlit] = ismember(psSeedIndexg,allSeedIndexg);
seedStFlit=gather(seedStFlit);
totalCnt=seedsWMax.weights(seedStFlit,:);



% [~,seedStFlit] = ismember(xSeedIndex,seedsWMax.seeds,'rows');
% totalCnt=seedsWMax.weights(seedStFlit,:);




sitesFlag=repelem((1:length(seedsCounts)).', seedsCounts,1);


seedsSites=vertcat(seedsSites{:});






nodes=seedsSites;
nodes=nodes.';
nodes=nodes(:);
    if length(nodes)>1e6
        warning('come on lad, input data is too big:%d in seedEnrichsCount \n',length(nodes));
    end


iSit=getZNICIM(nodes,  ones(W, 1));
iSit=sum(reshape(iSit, W, []));
iSit=iSit.';
iSit=iSit>0;



sitesFlagZp=sitesFlag(iSit);

pWeights=pWeights(iSit);

[sitesU, iz, jz]=uniqueSorted(sitesFlagZp);


nuMers=length(iz);
zoopsCntT=zeros(nuMers, 1);
if ~isempty(jz)
    pCountst = accumarray(jz,pWeights);
    zoopsCntT(1:length(pCountst))=pCountst;
end

zoopsCnt=zeros(numSeeds, 1);
zoopsCnt(sitesU)=zoopsCntT;

