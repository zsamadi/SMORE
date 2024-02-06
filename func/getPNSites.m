function seedsWmax=getPNSites(seqData,seedsWmax, options)

% Testing refined motifs on the holdout data
% input arguments:
%   - seedsRefined: struct for sorted input WMers
%       - seedsRefined.seeds   :  Refined NREF consensus motifs
%       - seedsRefined.pvalues : significance obtained by optimum
%       thresholding
%       - seedsRefined.PWMSCell: PWMS motif to be used for enrichment test
%       - seedsRefined.thresholdOptimum  : optimum thresholds
%   - seqData:
%       - seqData.pSeq  : primary sequence
%       - seqData.nSeq  : control sequence
%       - seqData.pHSeq : primary holdoput sequence
%       - seqData.nHSeq : control holdoput sequence
%       - seqData.PWM0  : backgroung model
%       - seqData.lens  : data lengths [size(posSeq90, 1),size(negSeq90, 1),size(pHoldSeq, 1),size(nHoldSeq, 1)];
%  Output arguments:
%   - outMotif: struct for sorted input WMers
%       - outMotif.consensusSeed   :  output consensus motif
%       - outMotif.testPvalue   : significance obtained from hold-out
%       - outMotif.trainPvalue  : significance obtained from training
%       - outMotif.PWMSE           : out PWMS motif to be erased
%       - outMotif.scoreThr        : Score thresholds



% G=options.G;


options.isPos=true;
options.isPN=false;



unqMerP=getPSNmers(seqData, options);

seedsWmax.posPWMS=unqMerP.posPWMS;

seedsWmax.pSSeeds=unqMerP.seeds;
seedsWmax.pSiteCell=unqMerP.siteCell;
seedsWmax.pWeightCell=unqMerP.weightCell;



options.isPos=false;
% totalScale(1)=sum(unqMersP.counts(:, 1));

NWeightCellC=cell(options.gTrainNum, 1);
NSiteCellC=cell(options.gTrainNum, 1);
nSSeedsCell=cell(options.gTrainNum, 1);
negPWMSC=cell(options.gTrainNum, 1);

cNTypesMat=seqData.cNTypes;




for iter=1:options.gTrainNum
    seqData.cNTypes=cNTypesMat(:, iter);

    unqMerN=getPSNmers(seqData, options);


    NWeightCellC{iter}=unqMerN.weightCell;
    NSiteCellC{iter}=unqMerN.siteCell;
    nSSeedsCell{iter}=unqMerN.seeds;

    negPWMSC{iter}=unqMerN.posPWMS;


    

end

seedsWmax.negPWMS=negPWMSC;
seedsWmax.nSSeeds=nSSeedsCell;
seedsWmax.nSiteCell=NSiteCellC;
seedsWmax.nWeightCell=NWeightCellC;















