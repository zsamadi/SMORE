function [seqFiltered, erasedSites, erasedNodes, nodeErzFlagO]=seqFilterNew(seq,cTypes, filterSpecs)

% Erasing PWMS motif obtained from enrichment step with score above threshold
% Input arguments:
%   - seq       : input sequence
%   - PWMS      : Refined PWM Score matrix
%   - threshold : score threshod for erasing
% Output arguments:
%  - seqFiltered:input sequence with motif erased

% 
% filterSpecs.back=seqData.back;
% filterSpecs.order=seqData.order;
iSeqFilteredOut=filterSpecs.iSeqFilteredOut;

PWMS=filterSpecs.PWMS;
% threshold=filterSpecs.scoreThr;
threshold=0;

if (filterSpecs.rvp)
    numSeqs=size(seq,1)/2;


else
    numSeqs=size(seq,1);
end





[numCells, W]=size(PWMS);


cnvOptions.numPSeqs=numSeqs;
cnvOptions.W=W;
% cnvOptions.rvp=filterSpecs.rvp;
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

% sites=XNmerso(:, end-1);

% merEyesOnt=merEyesOn;
% XNmersot=XNmerso;

% if (filterSpecs.rvp)
%     rvpFlag=XNmerso(:, 1)>XNmerso(:, W);
%     XNmerso(rvpFlag, 1:W)=XNmerso(rvpFlag, W:-1:1);
%     merEyesOnW=(reshape(merEyesOn, W, [])).';
%     merEyesOnW(rvpFlag, :)=merEyesOnW(rvpFlag, end:-1:1);
%     merEyesOnW=merEyesOnW.';
%     merEyesOn=merEyesOnW(:);
% 
% end




[XNmerUo, ~, jxo]=unique(XNmerso(:, 1:W), 'rows');


% PWMScoreLettero=ones(size(XNmerUo));

% eraseFlago=(PWMScoreo>=threshold);
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
% erasedNodes=erasedNodes;

if iSeqFilteredOut

    % [PWMScoreo, PWMScoreLettero]=scoreWords(XNmerUo, PWMS, filterSpecs);
     PWMScoreLettero=ones(size(XNmerUo));
    PWMScoreLettero=PWMScoreLettero(jxo, :);
    
    
    
    eraseFlagLettero=(PWMScoreLettero>0);
    eraseFlagLettero=eraseFlagLettero&eraseFlago;
    
    eraseFlagLetteron=eraseFlagLettero.';
    
    eraseFlagLetteron=eraseFlagLetteron(:);
end

%Also erase reverse path
if filterSpecs.rvp

    
%     eraseFlagoRVP=(PWMScoreoRVP>=threshold);
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

    % [PWMScoreoRVP, PWMScoreLetteroRVP]=scoreWords(XNmerUo, PWMS(:, end:-1:1), filterSpecs);
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
    
    % erasedSites=reshape(eraseFlagLetteront, size(seq, 2), []);
    % 
    % erasedSites=erasedSites.';
    eraseLetters=eraseFlagLetteront>0;
    
    seqFilteredERS=seqH;
    
    % erasedNodes=seqFilteredERS(eraseLetters);
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


% seqNum=repelem((1:length(seqLens)).', seqLens-W+1);
% 
% erasedSites=accumarray(seqNum, eraseFlago|eraseFlagoRVP);


