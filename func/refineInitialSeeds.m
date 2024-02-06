
function seedsRefined=refineInitialSeeds(seedsNEVAL,seedsWmax,seqData, seqDataSpecs, options)

% from initial NEVAL motifs, chooses NREF motifs by counting all matches
% and approximate matches in primary and control data to the motifs and computing their pvalue
% input arguments:
%   - seedsInfo
%       - seedsInfo.seeds   : sorted seeds by pvalue
%       - seedsInfo.pvalues : sorted pvalues
%   - XnMerZoops: Cell containing zoops WMers of primary and control
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
%   - seedsRefined: struct for sorted input WMers
%      - seedsRefined.seeds: sorted NEVAL seeds, NREF of which to be further refined
%      - seedsRefined.pvalues: significance obtained by fisher exact test

prior=options.prior;
wMax=options.wMax;
NEVAL=options.NEVAL;


motifsAll=seedsNEVAL.seeds;
wSeedsAll=seedsNEVAL.wSeeds;
pvaluesAll=seedsNEVAL.pvalues;
maxWidthAll=seedsNEVAL.maxWidth;



wSeedsU=unique(wSeedsAll);
numW=length(wSeedsU);

wSeedsFlag=(wSeedsAll==wSeedsU.');

validEx=~any(motifsAll==options.numCells+1, 2);

motifsAllValid=motifsAll(validEx, :);
% motifsAllNotValid=motifsAll(~validEx, :);

pvaluesAllValid=pvaluesAll(validEx);
wSeedsAllValid=wSeedsAll(validEx);
maxWidthAllValid=maxWidthAll(validEx);
% wSeedsAllNotValid=wSeedsAll(~validEx);


wSeedsFlag=wSeedsFlag(validEx, :);


motifsNEVAL=zeros(numW*NEVAL,wMax);
wSeedsNEVAL=zeros(numW*NEVAL,1);
pvaluesNEVAL=zeros(numW*NEVAL,1);
maxWidthNEVAL=zeros(numW*NEVAL,1);


motifsNotNEVAL=cell(numW, 1);
wSeedsNotNEVAL=cell(numW, 1);
pvaluesNotNEVAL=cell(numW, 1);
maxWidthNotNEVAL=cell(numW, 1);
idx0=0;
for iW=1:numW

    motifsi=motifsAllValid(wSeedsFlag(:,iW), :);
    NEVAL0=min(NEVAL, size(motifsi, 1));
    motifsNEVAL(idx0+1:idx0+NEVAL0, :)=motifsi(1:NEVAL0, :);
    motifsNotNEVAL{iW}=motifsi(NEVAL+1:end, :);

    wSeedsAlli=wSeedsAllValid(wSeedsFlag(:,iW));
    wSeedsNEVAL(idx0+1:idx0+NEVAL0)=wSeedsAlli(1:NEVAL0);
    wSeedsNotNEVAL{iW}=wSeedsAlli(NEVAL0+1:end);

    pvaluesAllValidi=pvaluesAllValid(wSeedsFlag(:,iW));
    pvaluesNEVAL(idx0+1:idx0+NEVAL0)=pvaluesAllValidi(1:NEVAL0);
    pvaluesNotNEVAL{iW}=pvaluesAllValidi(NEVAL0+1:end);


    maxWidthAllValidi=maxWidthAllValid(wSeedsFlag(:,iW));
    maxWidthNEVAL(idx0+1:idx0+NEVAL0)=maxWidthAllValidi(1:NEVAL0);
    maxWidthNotNEVAL{iW}=maxWidthAllValidi(NEVAL0+1:end);

    idx0=idx0+NEVAL0;

end

NEVALSeeds=idx0;

pvaluesNEVAL=pvaluesNEVAL(1:NEVALSeeds);
motifsNEVAL=motifsNEVAL(1:NEVALSeeds, :);
wSeedsNEVAL=wSeedsNEVAL(1:NEVALSeeds);
maxWidthNEVAL=maxWidthNEVAL(1:NEVALSeeds);



[~, pidx]=sort(pvaluesNEVAL, 'ascend');
% pvaluesNEVAL=pvalmaxWidthNEVAL(:, 1);
motifsNEVAL=motifsNEVAL(pidx, :);
wSeedsNEVAL=wSeedsNEVAL(pidx);
maxWidthNEVAL=maxWidthNEVAL(pidx);


% motifsNotNEVAL=cell2mat(motifsNotNEVAL);
% motifsAll=[motifsNEVAL;motifsNotNEVAL;motifsAllNotValid];

% wSeedsNotNEVAL=cell2mat(wSeedsNotNEVAL);
% wSeedsAll=[wSeedsNEVAL;wSeedsNotNEVAL;wSeedsAllNotValid];





PWM0=seqDataSpecs.back{1};
priorBg=prior*PWM0/(1+prior);
PWM0L=log2(PWM0);



pvaluesLogEnrchd=zeros(NEVALSeeds, 1);

options.nRefIter=1;
options.isOnlyPScore=true;


for imt=1:NEVALSeeds

%     wSeed=wSeedsNEVAL(imt);

    motifSeq=motifsNEVAL(imt, :);
%         if all(motifSeq==[3, 5, 9, 10])
%             check=1;
%         end


    % initialize PWM1 with motif

    PWM1=(motifSeq(:)==(1:length(PWM0)))/(1+prior)+priorBg.';
    PWM1L=log2(PWM1.');
%     PWMS(:,1:wSeed)=(PWM1L(1:wSeed, :)-PWM0L).';
    if options.mkvOrder>0
        PWMS=PWM1L;
    else
        PWMS=(PWM1L-PWM0L);
    end
    options.seedNumber=[imt,NEVALSeeds];

     printString=sprintf('(%d,%d)', imt, NEVALSeeds);

     fprintf(printString);




    seedEnriched=nestedSeedEnrichment(PWMS, seedsWmax,seqData, seqDataSpecs, options);
    pvaluesLogEnrchd(imt)=seedEnriched.pvalue;

    delString=repelem(sprintf('\b'), 1, length(printString));

    fprintf(delString)



end




[pvalueNREFSort, sindex]=sort([pvaluesLogEnrchd, -wSeedsNEVAL],1);
pvalueNREFSort=pvalueNREFSort(:,1);
motifsNEVALNREF=motifsNEVAL(sindex, :);
wSeedsNEVAL=wSeedsNEVAL(sindex);
maxWidthNEVAL=maxWidthNEVAL(sindex);

[motifsNEVALNREF, iu]=unique(motifsNEVALNREF, 'stable', 'rows');
pvalueNREFSort=pvalueNREFSort(iu);
wSeedsNEVAL=wSeedsNEVAL(iu);
maxWidthNEVAL=maxWidthNEVAL(iu);

seedsRefined.seeds=motifsNEVALNREF;
seedsRefined.pvalues=pvalueNREFSort;
seedsRefined.wSeeds=wSeedsNEVAL;
seedsRefined.maxWidth=maxWidthNEVAL;
seedsRefined.seedsChar=[char(motifsNEVALNREF+64), num2str([pvalueNREFSort,maxWidthNEVAL, wSeedsNEVAL])];