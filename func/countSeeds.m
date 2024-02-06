function [unqMersCell, unqMersFlagsCell]=countSeeds(seqData, unqMersPIn,options)


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

MinSeedWidth=options.MinSeedWidth;
wMax=options.wMax;
numCells=options.numCells;
numNodes=length(seqData.cPTypes);

cTypesEx=[seqData.cPTypes; numCells+1;0;seqData.cNTypes];
cnvOptions.cTypesInit=[options.cPTypesInit; numCells+1;0;options.cNTypesInit];


% Sorts input sequence matrix in NMers and Zoops NMers
lenExt=wMax-MinSeedWidth;


pSeq=seqData.pSeq;
nSeq=seqData.nSeq;

numPSeqs=length(pSeq);

pExt=repelem(numNodes+1,1, lenExt);

pSeqE=cell(seqData.lens(1),1);
pSeqLens=zeros(seqData.lens(1), 1);
for iPs=1:numPSeqs
    pSeqEi=[pSeq{iPs},pExt];
    pSeqE{iPs}=pSeqEi;
    pSeqLens(iPs)=length(pSeqEi);
end

pSeqESt=horzcat(pSeqE{:});


nSeqE=cell(seqData.lens(2),1);
nSeqLens=zeros(seqData.lens(2), 1);
numNSeqs=length(nSeq);
nExt=repelem(-1,1, lenExt);
for iNs=1:numNSeqs
    nSeqEi=[nSeq{iNs},nExt]+numNodes+2;
    nSeqE{iNs}=nSeqEi;
    nSeqLens(iNs)=length(nSeqEi);

end



    



nSeqESt=horzcat(nSeqE{:});





pnseqEx=[pSeqESt,nSeqESt];
pnSeqLens=[pSeqLens;nSeqLens];

if options.rvp
    cnvOptions.numPSeqs=numPSeqs/2;
else
    cnvOptions.numPSeqs=numPSeqs;
end
cnvOptions.W=wMax;
cnvOptions.rvp=options.rvp;

cnvOptions.allMers=false;
cnvOptions.pnSeq=true;
cnvOptions.numCells=options.numCells;
cnvOptions.erasePhase=false;

cnvOptions.isOLess=options.isOLess;
cnvOptions.isOKSingleNC=options.isOKSingleNC;



nMers=convertNMer(pnseqEx,pnSeqLens,cTypesEx, cnvOptions);




% if any(pnSeqLens>wMax+1)
%     aa=nMers(:,wMax+1:end-1);
%     [~, iU]=unique(aa, 'rows');
%     nMers=nMers(iU, :);
% end

nMers(:, end)=nMers(:, end)-lenExt;



numW=wMax-MinSeedWidth+1;
unqMersCell=cell(numW, 1);
unqMersFlagsCell=cell(numW, 1);

options.back=seqData.back;



% cutEnd=wMax-MIN_SEED_WIDTH;
pnWeights=[seqData.posWeight;seqData.negWeight];
repNum=size(nMers, 1)/length(pnWeights);

if repNum>1
    pnWeights=repelem(repNum, repNum, 1);
end


options.pnWeights=pnWeights;


wSeedv=(MinSeedWidth:wMax);
unqSeedE=[];

for wseedi=numW:-1:1
    wseed=wSeedv(wseedi);
    % [unqMers,unqSeedE, unqMersFlagsOut] =countWSeeds(nMers, wseed,unqSeedE,unqMersFlagsIn{wseedi}, options);
    [unqMers,unqSeedE] =countWSeeds(nMers, wseed,unqSeedE,unqMersPIn{wseedi}, options);

    unqMersCell{wseedi}=unqMers;
    % unqMersFlagsCell{wseedi}=unqMersFlagsOut;

end







