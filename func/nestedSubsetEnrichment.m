
function motifsEnriched=nestedSubsetEnrichment(seedsNREF, seedsWmax, seqData,seqDataSpecs, options)

NREF=options.NREF;
prior=options.REFPrior;

PWM0=seqDataSpecs.back{1};
NREF=min(NREF, size(seedsNREF.seeds, 1));
seedsNREF=seedsNREF.seeds(1:NREF, :);
priorBg=prior*PWM0/(1+prior);
PWM0L=log2(PWM0);
pvalueOptimum=ones(NREF,1);
thrOptVec=zeros(NREF,1);
PWMSCell=cell(NREF, 1);
PWMSOutCell=cell(NREF, 1);
seedsToPWMCell=cell(NREF, 1);
seedsPvalueCell=cell(NREF, 1);
consensusSeedMat=seedsNREF;

mkvOrder=options.mkvOrder;

for imotifSeq=1:NREF
    seedSeq =seedsNREF(imotifSeq, :);

    PWM1=(seedSeq(:)==(1:length(PWM0)))/(1+prior)+priorBg.';
    PWM1L=log2(PWM1.');

    if mkvOrder>0
        PWMS=PWM1L;
    else
        PWMS=(PWM1L-PWM0L);
    end

    options.seedNumber=[imotifSeq,NREF];
    printString=sprintf('(%d,%d)', imotifSeq, NREF);
    fprintf(printString);
    seedEnriched=nestedSeedEnrichment(PWMS, seedsWmax,seqData, seqDataSpecs, options);
    PWMSCell{imotifSeq}=seedEnriched.PWMS;
    PWMSOutCell{imotifSeq}=seedEnriched.PWMOut;
    thrOptVec(imotifSeq)=seedEnriched.thresholdOptimum;
    pvalueOptimum(imotifSeq)=seedEnriched.pvalue;
    consensusSeedMat(imotifSeq, :)=seedEnriched.seed;
    seedsToPWMCell{imotifSeq}=seedEnriched.seedsToPWM;
    seedsPvalueCell{imotifSeq}=seedEnriched.seedsPvalue;

    delString=repelem(sprintf('\b'), 1, length(printString));

    fprintf(delString)


end

[pvalueOptimum,idx]=sort(pvalueOptimum, 'ascend');
PWMSCell=PWMSCell(idx);
PWMSOutCell=PWMSOutCell(idx);
thrOptVec=thrOptVec(idx);
seedsNREF=seedsNREF(idx, :);
consensusSeedMat=consensusSeedMat(idx, :);
seedsToPWMCell=seedsToPWMCell(idx);
seedsPvalueCell=seedsPvalueCell(idx);

motifsEnriched.seeds=consensusSeedMat;
motifsEnriched.pvalues=pvalueOptimum;
motifsEnriched.PWMSCell=PWMSCell;
motifsEnriched.PWMSOutCell=PWMSOutCell;
motifsEnriched.thresholdOptimum=thrOptVec;
motifsEnriched.seedsChar=[char(seedsNREF+64), num2str([pvalueOptimum, thrOptVec])];
motifsEnriched.isPal=false;
motifsEnriched.seedsToPWM=seedsToPWMCell;
motifsEnriched.seedsPvalue=seedsPvalueCell;

if (options.isPal) && ~isempty(PWMSCell)
    [palModel, edOptimum]=getPalMode(PWMSCell{1}, PWM0);
    if seqDataSpecs.mkvOrder>0
        PWMSPal=log2(palModel);
    else

        PWMSPal=log2(palModel)-PWM0L;
    end

    options.isPal=true;
    motifPal=nestedSeedEnrichment(PWMSPal, seedsWmax, seqData,seqDataSpecs, options);
    logPvalueRatio=motifPal.pvalue/pvalueOptimum(1);
    if edOptimum<options.maxPalED && logPvalueRatio>=options.minPalRatio

        motifsEnriched.seeds(1, :)=motifPal.seed;
        motifsEnriched.pvalues(1)=motifPal.pvalue;
        motifsEnriched.PWMSCell{1}=motifPal.PWMS;
        motifsEnriched.PWMSOutCell{1}=motifPal.PWMOut;
        motifsEnriched.thresholdOptimum(1)=motifPal.thresholdOptimum;
        motifsEnriched.isPal=true;
        motifsEnriched.seedsToPWM{1}=motifPal.seedsToPWM;
        motifsEnriched.seedsPvalue{1}=motifPal.seedsPvalue;
    end
end



