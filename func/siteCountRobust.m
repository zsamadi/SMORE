function [posCntN, seedsInPWM]=siteCountRobust(unqMers, seedsToPWM)


W=size(unqMers.seeds, 2);
allWeights=unqMers.weightCell;
pnSites=unqMers.siteCell;
numSeeds=size(seedsToPWM, 1);
pnSitesIn=cell(numSeeds, 1);
if ~isempty(unqMers.seeds)
    [Lia, indx] = ismember(unqMers.seeds,seedsToPWM,'rows');

    pnSitesIn(indx(Lia))=pnSites(Lia);
    allWeightsIn=cell(numSeeds, 1);

    allWeightsIn(indx(Lia))=allWeights(Lia);

    seedsInPWMCounts=zeros(length(allWeightsIn), 1);
    for iAW=1:length(allWeightsIn)
        if ~isempty(allWeightsIn{iAW})
            seedsInPWMCounts(iAW)=sum(allWeightsIn{iAW}, 1);
        end
    end


    seedsInPWM=[seedsToPWM, seedsInPWMCounts(:, 1)];


    pnSitesIn=vertcat(pnSitesIn{:});
    allWeights=vertcat(allWeights{:});


    cNodes=pnSitesIn.';
    cNodes=cNodes(:);
    if ~isempty(cNodes)
        unqMersFlags=getZNICIM(cNodes,  ones(W, 1));
        unqMersFlags=reshape(unqMersFlags, W, []);
        unqMersFlags=sum(unqMersFlags);
        unqMersFlags=unqMersFlags.';

        allWeights=allWeights(unqMersFlags>0);

        posCntN=sum(allWeights);
    else
        posCntN=0;
    end

else
    posCntN=0;
    seedsInPWM=[seedsToPWM, zeros(numSeeds, 1)];
end

