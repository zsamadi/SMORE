function [posSeq, pHoldSeq, cPTypes, pWeight, pHdWeight]=generateSeqs(pData,pWeight, dPSections, cPTypes,options)
% from input sequence, generates control sequence by shuffling, and
% separates 10% of the data as hold out data

% input arguments:
% - dstream: input sequence
% % Output arguments:
% - posSeq90: output primary sequence(90% of input sequence)
% - negSeq90: output control sequence(90% of shuffled input sequence)
% - pHoldSeq: output primary hold out sequence
% - nHoldSeq: output control hold out sequence



hFrac=options.hFrac;
order=options.mkvOrder;

rvpath=options.rvp;

numPSeqs=length(pData);



if options.numGRs>1
     pDataM=vertcat(pData{:});
     nodeINC=options.numNodes*(0:options.numGRs-1).';
     nodeINC=repelem(nodeINC,numPSeqs,1);
     pDataM50=repmat(pDataM, options.numGRs,1);
     pDataM50=pDataM50+nodeINC;
     pData=num2cell(pDataM50,2);
     pWeight=repmat(pWeight, options.numGRs, 1);
end
numPSeqs=length(pData);

if (rvpath)
    dstreamInv=cell(numPSeqs, 1);
    for iSp=1:numPSeqs
        dstreamInv{iSp}=pData{iSp}(end:-1:1);
    end

end

pHoldIDX=false(numPSeqs,1);

if options.isSectHold && hFrac>0
    pDSections=dPSections;
   pHoldIDX=(pDSections==1);

else
    numPHSeq=floor(numPSeqs*hFrac);
    randIDX=randperm(numPSeqs);
    pHoldIDX(randIDX(1:numPHSeq))=true;

end





pHoldSeq=pData(pHoldIDX);
pHdWeight=pWeight(pHoldIDX);


pHoldSeqAll=vertcat(pHoldSeq{:});
pHoldSeqAll=pHoldSeqAll(:);

if ~options.isSectHold
    cPTypes(pHoldSeqAll)=0;
end



if (rvpath)

    if options.isHExclsv
         posSeq=vertcat(pData(~pHoldIDX), dstreamInv(~pHoldIDX));
         pWeight=repmat(pWeight(~pHoldIDX), 2, 1);
    else
         posSeq=vertcat(pData, dstreamInv);
         pWeight=repmat(pWeight, 2, 1);

    end

    pHoldSeq=vertcat(pHoldSeq, dstreamInv(pHoldIDX));
    pHdWeight=repmat(pHdWeight, 2, 1);

end


if isempty(pHoldSeq)
    pHoldSeq=posSeq;
    pHdWeight=pWeight;
end



check=1;





