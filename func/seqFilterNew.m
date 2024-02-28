function [seqFiltered, erasedSites, erasedNodes, nodeErzFlagO]=seqFilterNew(seq,cTypes, filterSpecs)

iSeqFilteredOut=filterSpecs.iSeqFilteredOut;

PWMS=filterSpecs.PWMS;
threshold=0;

if (filterSpecs.rvp)
    numSeqs=size(seq,1)/2;
else
    numSeqs=size(seq,1);
end

[numCells, W]=size(PWMS);
cnvOptions.numPSeqs=numSeqs;
cnvOptions.W=W;
cnvOptions.rvp=false;
cnvOptions.allMers=true;
cnvOptions.pnSeq=false;
cnvOptions.numCells=numCells;
cnvOptions.erasePhase=true;
cnvOptions.isOLess=filterSpecs.isOLess;
cnvOptions.cTypesInit=filterSpecs.cPTypesInit;

cnvOptions.isOKSingleNC=(filterSpecs.isOKSingleNC && filterSpecs.isEraseNodes);
cnvOptions.isNotNMer=filterSpecs.isNotNMer;

if filterSpecs.isNotNMer

    seqH=horzcat(seq{1:numSeqs});

    seqLens=zeros(numSeqs, 1);
    for iSl=1:numSeqs
        seqLens(iSl)=length(seq{iSl});
    end

    seqLens=zeros(numSeqs, 1);
    for iPs=1:numSeqs
        seqLens(iPs)=length(seq{iPs});
    end
else
    seqH=seq(1:numSeqs, :);
    seqLens=repelem(W,numSeqs, 1);

end

[XNmerso , merEyesOn]=convertNMer(seqH,seqLens,cTypes, cnvOptions);
[XNmerUo, ~, jxo]=unique(XNmerso(:, 1:W), 'rows');


[seedPWMScoreo, seedPWMScoreLettero]=scoreWords(filterSpecs.seedsToPWM, PWMS, filterSpecs);
seedPWMScoreLettero=seedPWMScoreLettero(seedPWMScoreo>=threshold, :);
seedLetterErzFlag=seedPWMScoreLettero>0;
seedsToPWM=filterSpecs.seedsToPWM(seedPWMScoreo>=threshold, :);

[eraseFlago, ila]=ismember(XNmerUo, seedsToPWM, 'rows');

nodeErzFlagO=false(size(XNmerUo));
nodeErzFlagO(eraseFlago, :)=seedLetterErzFlag(ila(ila>0), :);
nodeErzFlagO=nodeErzFlagO(jxo, :);

eraseFlago=eraseFlago(jxo);
nodeErzFlagO=nodeErzFlagO(eraseFlago, :);

erasedNodes=XNmerso(eraseFlago, :);

erasedNodes=erasedNodes(:, W+1:2*W);

if iSeqFilteredOut

    PWMScoreLettero=ones(size(XNmerUo));
    PWMScoreLettero=PWMScoreLettero(jxo, :);



    eraseFlagLettero=(PWMScoreLettero>0);
    eraseFlagLettero=eraseFlagLettero&eraseFlago;

    eraseFlagLetteron=eraseFlagLettero.';

    eraseFlagLetteron=eraseFlagLetteron(:);
end

if filterSpecs.rvp


    [eraseFlagoRVP, ilrvp]=ismember(XNmerUo, seedsToPWM(:, end:-1:1), 'rows');

    rvSeedLetterErzFlag=seedLetterErzFlag(:, end:-1:1);

    rvSeedLetterErzFlagO=false(size(XNmerUo));
    rvSeedLetterErzFlagO(eraseFlagoRVP, :)=rvSeedLetterErzFlag(ilrvp(ilrvp>0), :);
    rvSeedLetterErzFlagO=rvSeedLetterErzFlagO(jxo, :);



    eraseFlagoRVP=eraseFlagoRVP(jxo);
    rvSeedLetterErzFlagO=rvSeedLetterErzFlagO(eraseFlagoRVP, :);
    nodeErzFlagO=[nodeErzFlagO;rvSeedLetterErzFlagO];

    erasedNodesRv=XNmerso(eraseFlagoRVP, :);
    erasedNodesRv=erasedNodesRv(:, 2*W:-1:W+1);
    erasedNodes=[erasedNodes;erasedNodesRv];

    if iSeqFilteredOut

        PWMScoreLetteroRVP=ones(size(XNmerUo));
        PWMScoreLetteroRVP=PWMScoreLetteroRVP(jxo, :);
        eraseFlagLetteroRVP=(PWMScoreLetteroRVP>0);
        eraseFlagLetteroRVP=eraseFlagLetteroRVP&eraseFlagoRVP;
        eraseFlagLetteronRVP=eraseFlagLetteroRVP.';

        eraseFlagLetteronRVP=eraseFlagLetteronRVP(:);
        eraseFlagLetteron=eraseFlagLetteron|eraseFlagLetteronRVP;
    end

end


if iSeqFilteredOut

    seqH=seqH(:, 1:W);
    seqH=seqH.';
    seqH=seqH(:);

    eraseFlagLetteront=zeros(numel(seqH), 1);
    temp=accumarray(merEyesOn, eraseFlagLetteron);

    eraseFlagLetteront(1:length(temp))=temp;

    eraseLetters=eraseFlagLetteront>0;

    seqFilteredERS=seqH;

    seqFilteredERS(eraseLetters)=length(cTypes);
    seqFiltered=reshape(seqFilteredERS, W,[]);
    seqFiltered=seqFiltered.';

    if (filterSpecs.rvp)
        seqFilteredRVP=reshape(seqFilteredERS(end:-1:1),W, []);
        seqFilteredRVP=seqFilteredRVP.';
        seqFilteredRVP=seqFilteredRVP(end:-1:1, :);
        seqFiltered=[seqFiltered;seqFilteredRVP];
        seqFiltered=[seqFiltered, seq(:, end)];
    end
else
    seqFiltered=[];
end

erasedSites=eraseFlago|eraseFlagoRVP;


