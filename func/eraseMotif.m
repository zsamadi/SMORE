function [seqData, sitesErased, isEraseOccured]=eraseMotif(outMotif,  seqData, options)

% Erasing final motif obtained from enrichment step
% Input arguments:
%  - outMotif:
%       - outMotif.PWMSE        : Final refined PWM score matrix
%       - outMotif.consensusSeed: Final consensus Seed
%       - outMotif.scoreThr     : score threshold used for erasing
%       - outMotif.enrichPvalue : Enrichment pvalue of the output motif
%   - seqData: 
%       - seqData.pSeq  : primary sequence
%       - seqData.nSeq  : control sequence
%       - seqData.pHSeq : primary holdoput sequence
%       - seqData.nHSeq : control holdoput sequence
%       - seqData.PWM0  : backgroung model
%       - seqData.lens  : data lengths [size(posSeq90, 1),size(negSeq90, 1),size(pHoldSeq, 1),size(nHoldSeq, 1)];
% Output arguments:
%   - seqData: output seqData structure with input motif erased



% Filter the input sequence with threshold score, letters only erased if
% they're part of a motif and their score is greater than zero

filterSpecs=options;


filterSpecs.PWMS=outMotif.PWMSE;
filterSpecs.scoreThr=outMotif.scoreThr;
filterSpecs.back=seqData.back;
filterSpecs.seedsToPWM=outMotif.seedsToPWM;

filterSpecs.iSeqFilteredOut=~options.isEraseNodes;

cTypesEx=[seqData.cPTypes(:); 0;0];

[pSeqFiltered, pSeqErasedSites, pSeqErasedNodes, pNodeErzFlag]=seqFilterNew(seqData.pSeq,cTypesEx, filterSpecs);

filterSpecs.iSeqFilteredOut=false;

if options.isErzHOut
    cTypesEx=[seqData.cPHTypes(:); 0;0];
    [pHSeqFiltered,pHSeqErasedSites, pHSeqErasedNodes, pHNodeErzFlag]=seqFilterNew(seqData.pSeq,cTypesEx, filterSpecs);
else
   pHSeqFiltered= seqData.pSeq;
   pHSeqErasedSites=0;
   pHSeqErasedNodes=[];
   pHNodeErzFlag=[]>0;
end

if options.isErzNegative
    filterSpecs.cTypesInit=options.cNTypesInit;
    cTypesEx=[seqData.cNTypes(:); 0;0];
    [nSeqFiltered, nSeqErasedSites, nSeqErasedNodes, nNodeErzFlag]=seqFilterNew(seqData.pSeq,cTypesEx, filterSpecs);

    if options.isErzHOut
        cTypesEx=[seqData.cNHTypes(:); 0;0];
        [nHSeqFiltered,nHSeqErasedSites, nHSeqErasedNodes, nHNodeErzFlag]=seqFilterNew(seqData.pSeq,cTypesEx, filterSpecs);
    else
       nHSeqFiltered= seqData.nHSeq;
       nHSeqErasedSites=0;
       nHSeqErasedNodes=[];
       nHNodeErzFlag=[]>0;
    end    
else
    nSeqFiltered=seqData.pSeq;
    nSeqErasedSites=0;
    nHSeqFiltered=seqData.pSeq;
    nHSeqErasedSites=0;
    nSeqErasedNodes=[];
    nHSeqErasedNodes=[];
    nNodeErzFlag=[]>0;
    nHNodeErzFlag=[]>0;
end


% seqDataCopy.pSeq=seqData.pSeq;
% seqDataCopy.nSeq=seqData.nSeq;
% seqDataCopy.pHSeq=seqData.pHSeq;
% seqDataCopy.nHSeq=seqData.nHSeq;


if ~options.isEraseNodes
    seqData.pSeq=pSeqFiltered;
    % seqData.nSeq=nSeqFiltered;
    % seqData.pHSeq=pHSeqFiltered;
    % seqData.nHSeq=nHSeqFiltered;
end

sitesErased.pSeq=pSeqErasedSites;
sitesErased.nSeq=nSeqErasedSites;
sitesErased.pHSeq=pHSeqErasedSites;
sitesErased.nHSeq=nHSeqErasedSites;


sitesErased.pNodes=pSeqErasedNodes;
sitesErased.nNodes=nSeqErasedNodes;
sitesErased.pHNodes=pHSeqErasedNodes;
sitesErased.nHNodes=nHSeqErasedNodes;

pNodes=pSeqErasedNodes;
pHNodes=pHSeqErasedNodes;


ppHNode=[pNodes;pHNodes];

ppHNode=ppHNode(:);
% pNodes=pNodes(:);
pHNodes=pHNodes(:);
% ppHNode=unique(ppHNode);
% seqData.cPTypes(ppHNode)=0;
        % seqData.cPHTypes(ppHNode)=0;

nnHNode=[nSeqErasedNodes;nHSeqErasedNodes];
nnHNode=nnHNode(:);


    pNodeFlag=pNodeErzFlag;
    pHNodeFlag=pHNodeErzFlag;

    ppHNodeFlag=[pNodeFlag;pHNodeFlag];
    ppHNodeFlag=ppHNodeFlag(:);
    ppHNode=ppHNode(ppHNodeFlag);

    pHNodes=pHNodes(pHNodeFlag(:));




    nnHNodeFlag=[nNodeErzFlag;nHNodeErzFlag];
    nnHNodeFlag=nnHNodeFlag(:);
    nnHNode=nnHNode(nnHNodeFlag);

    ppHNode=unique(ppHNode);
    pHNodes=unique(pHNodes);



if options.isErzFixNodes
    shConfig=options.shConfig;
    seqData.fixedNodes=[seqData.fixedNodes;ppHNode];
    shConfig.fixedNodes=seqData.fixedNodes;
    seqData.cNTypes=getShuffleTypes(seqData.cPTypes, shConfig);
else

    if options.isEraseNodes
        seqData.cPTypes(ppHNode)=0; 
        seqData.cPHTypes(pHNodes)=0; 
    
    
        nnHNode=unique(nnHNode);
        seqData.cNTypes(nnHNode)=0;
    end
end

isEraseOccured=any(pSeqErasedSites);

