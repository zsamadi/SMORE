function [seqData, sitesErased, isEraseOccured]=eraseMotif(outMotif,  seqData, options)


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
    [~,pHSeqErasedSites, pHSeqErasedNodes, pHNodeErzFlag]=seqFilterNew(seqData.pSeq,cTypesEx, filterSpecs);
else
    pHSeqErasedSites=0;
    pHSeqErasedNodes=[];
    pHNodeErzFlag=[]>0;
end

if options.isErzNegative
    filterSpecs.cTypesInit=options.cNTypesInit;
    cTypesEx=[seqData.cNTypes(:); 0;0];
    [~, nSeqErasedSites, nSeqErasedNodes, nNodeErzFlag]=seqFilterNew(seqData.pSeq,cTypesEx, filterSpecs);

    if options.isErzHOut
        cTypesEx=[seqData.cNHTypes(:); 0;0];
        [~,nHSeqErasedSites, nHSeqErasedNodes, nHNodeErzFlag]=seqFilterNew(seqData.pSeq,cTypesEx, filterSpecs);
    else

        nHSeqErasedSites=0;
        nHSeqErasedNodes=[];
        nHNodeErzFlag=[]>0;
    end
else
    nSeqErasedSites=0;
    nHSeqErasedSites=0;
    nSeqErasedNodes=[];
    nHSeqErasedNodes=[];
    nNodeErzFlag=[]>0;
    nHNodeErzFlag=[]>0;
end



if ~options.isEraseNodes
    seqData.pSeq=pSeqFiltered;
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
pHNodes=pHNodes(:);


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

