function [unqMers, seqData]=getNMer(seqData, options)


% Counts all WMers in the input primary and control data
% input arguments:
%   - seqData:
%       - seqData.pSeq  : primary sequence
%       - seqData.nSeq  : control sequence
%       - seqData.pHSeq : primary holdoput sequence
%       - seqData.nHSeq : control holdoput sequence
%       - seqData.PWM0  : backgroung model
%       - seqData.lens  : data lengths [size(posSeq90, 1),size(negSeq90, 1),size(pHoldSeq, 1),size(nHoldSeq, 1)];
%   - options
%       - options.NEVAL     : Number of motifs for initial evaluation
%       - options.NREF      : Number of motifs for refinement
%       - options.prior     : Dirichlet prior for initial evaluation
%       - options.REFPrior  : Dirichlet prior for refinement step
%  Output arguments:
%   - infoCell
%     - uniMers     : all seeds of length W, counts of Wmer seeds in primary and control data
%     - primCntsF   : counts of Wmer seeds in primary data
%     - ctrlCntsF   : counts of Wmer seeds in control data
%   - wmerZoops     : structure of cells containing zoops WMers of primary and control
%       - wmerZoops.pos=posWmerZoopsCell;
%       - wmerZoops.neg=negWmerZoopsCell;
% wMin=options.wMin;
wMax=options.wMax;
numCells=options.numCells;

isUHold=options.isUHold;

numPSeqs=length(seqData.pSeq);

cnvOptions.W=wMax;
cnvOptions.allMers=false;
cnvOptions.isOLess=options.isOLess;
cnvOptions.isOKSingleNC=options.isOKSingleNC;








if options.isPos


    if isUHold
        cHTypesEx=[seqData.cPHTypes(:); numCells+1;0];
        cnvOptions.cTypesInit=[options.cPTypesInit; numCells+1;0];

    else
        cHTypesEx=[seqData.cPTypes(:); numCells+1;0];
        cnvOptions.cTypesInit=[options.cPTypesInit; numCells+1;0];


    end




else


    if isUHold
        cHTypesEx=[seqData.cNHTypes(:); numCells+1;0];
        pWeights=seqData.nHoldWeight;
        cnvOptions.cTypesInit=[options.cNTypesInit; numCells+1;0];


    else
        cHTypesEx=[seqData.cNTypes(:); numCells+1;0];
        cnvOptions.cTypesInit=[options.cNTypesInit; numCells+1;0];

    end


end
% Sorts input sequence matrix in NMers and Zoops NMers




    pSeqLens=zeros(numPSeqs, 1);
    for iPs=1:numPSeqs
        pSeqLens(iPs)=length(seqData.pSeq{iPs});
    end
    
    pSeqESt=horzcat(seqData.pSeq{:});

    cnvOptions.isNotNMer=true;

    nMers=convertNMer(pSeqESt,pSeqLens,cHTypesEx, cnvOptions);
    seqData.pSeq=nMers(:, wMax+1:end);











