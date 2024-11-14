function [seedEnriched, iterNumber]=nestedSeedEnrichment(PWMS, seedsWMax,seqData, seqDataSpecs, options)


refIter=1;
allPvalues=seedsWMax.pvalues;

allSeeds=seedsWMax.seeds;

addPos=options.addPos;

PWM0=seqDataSpecs.back{1};
PWM0L=log2(PWM0);

nPos=floor((options.numNodes));
nNeg=nPos;

bernoulli=options.bernoulli;

prior=options.REFPrior;
nRefIter=options.nRefIter;

thrOptimum=0;
pvalue=2;
[numCells,wMax]=size(PWMS);
cellCharsDouble=(1:numCells);

if seqDataSpecs.mkvOrder>0
    PWM1=2.^(PWMS);
else

    PWM1=2.^(PWMS+PWM0L);
end
PWM10=PWM1;
pvalue0=1;
D15=numCells.^(wMax-1:-1:0);

[~, seedsToPWM]=max(PWM1);
seedsPvalue0=-1;
appxSeedsLogPvalues=-1;
options.isPN=false;

while (refIter<=nRefIter)

    posPWMS=scoreWords(allSeeds,PWMS, seqDataSpecs);
    positiveFlags=posPWMS>0;


    appxPvalues=allPvalues(positiveFlags);
    possibleSeeds=allSeeds(positiveFlags, :);

    positivePWMS=posPWMS(positiveFlags);

    [appxPWMS, iN]=sort(positivePWMS, 'descend');

    appxSeeds=possibleSeeds(iN, :);

    appxPvalues=appxPvalues(iN);

    if isempty(appxSeeds)
        pvalue=1;
        refIter0=refIter;
        refIter=nRefIter+1;
    elseif ~options.isEnrich
        pvalue=appxPvalues(1);
        refIter0=1;
        refIter=nRefIter+1;  
    else
        PWMSTV=appxPWMS;


        if size(appxSeeds, 1)>1
            options.PWMS=PWMS;
            options.isPos=true;
            cPTypes=seqData.cPTypes;

            unqMerP=getPSNmers(seqData,cPTypes, options);

            seedsWMax.posPWMS=unqMerP.posPWMS;

            seedsWMax.pSeeds=unqMerP.seeds;
            seedsWMax.pSiteCell=unqMerP.siteCell;
            seedsWMax.pWeightCell=unqMerP.weightCell;
            if ~isempty(unqMerP.seeds)
            seedsWMax.pSeedsIndex=(sum(unqMerP.seeds.*D15,2));
            appxSeedIndex=(sum(appxSeeds.*D15,2));



            [zoopsCntAllP,totalCntAllP] =seedEnrichsCount(seedsWMax,appxSeedIndex);
            totalCntAllP=totalCntAllP(:, 1);
            options.isPos=false;
            cNTypesMat=options.cNTypesMat;

            totalCntAllNS=zeros(size(appxSeeds, 1),options.gTrainNum+addPos);

            zoopsCntAllNS=zeros(size(appxSeeds, 1),options.gTrainNum+addPos);

            zoopsCntAllNS(:,1)=zoopsCntAllP;
            totalCntAllNS(:,1)=totalCntAllP;
            for iter=1:options.gTrainNum

                cNTypesi=cNTypesMat(:,iter);

                unqMerP=getPSNmers(seqData,cNTypesi, options);

                if ~isempty(unqMerP.seeds)
                    seedsWMaxi=seedsWMax;

                    seedsWMaxi.posPWMS=unqMerP.posPWMS;

                    seedsWMaxi.pSeeds=unqMerP.seeds;

                    seedsWMaxi.pSeedsIndex=(sum(unqMerP.seeds.*D15,2));

                    seedsWMaxi.pSiteCell=unqMerP.siteCell;
                    seedsWMaxi.pWeightCell=unqMerP.weightCell;

                    [zoopsCntAllN,totalCntAllN] =seedEnrichsCount(seedsWMaxi,appxSeedIndex);
                    totalCntAllN=totalCntAllN(:, 2);
                    zoopsCntAllNS(:,iter+addPos)=zoopsCntAllN;
                    totalCntAllNS(:,iter+addPos)=totalCntAllN;
                end

            end


            zoopsCntAllNS2=zoopsCntAllNS;
            zoopsCntAllNS=sum(zoopsCntAllNS, 2);
            totalCntAllNS=sum(totalCntAllNS, 2);

            zoopsCntAll=[zoopsCntAllP, zoopsCntAllNS];
            totalCntAll=[totalCntAllP, totalCntAllNS];

            effectIndex=sum(zoopsCntAll, 2)>0;

            appxSeeds=appxSeeds(effectIndex, :);
            appxPvalues=appxPvalues(effectIndex);
            zoopsCntAll=zoopsCntAll(effectIndex, :);
            totalCntAll=totalCntAll(effectIndex, :);
            zoopsCntAllNS2=zoopsCntAllNS2(effectIndex, :);

            appxPWMS=appxPWMS(effectIndex);
            if options.indSeedMode
                randDisturb=0.0001*rand(size(appxPWMS));
                randDisturb=sort(randDisturb);
            else
                randDisturb=0;
            end

            appxPWMS=appxPWMS-randDisturb;
            appxPWMS(appxPWMS<0)=0;

            [PWMSTV,~, jPMSU]=unique(appxPWMS);

            seedCntsThr= [accumarray(jPMSU,zoopsCntAll(:,1)),accumarray(jPMSU,zoopsCntAll(:,2))] ;
            seedCntsThr=cumsum(seedCntsThr(end:-1:1, :), 1);

            seedCntsThr(:, 2)=seedCntsThr(:, 2)/(options.gTrainNum+addPos);

            zoopsCntAllNS2M=zeros(length(PWMSTV), options.gTrainNum+addPos);

            for iACM=1:options.gTrainNum+addPos
                zoopsCntAllNS2M(:, iACM)=accumarray(jPMSU, zoopsCntAllNS2(:, iACM));
            end

            zoopsCntAllNS2M=cumsum(zoopsCntAllNS2M(end:-1:1, :), 1);

            if options.isNormalTest
                seedCntsThrVar= sum(zoopsCntAllNS2M.^2, 2)/(options.gTrainNum+addPos);

                lgnMu=log(seedCntsThr(:, 2).^2./sqrt(seedCntsThrVar));
                lgnSigma=sqrt(log(seedCntsThrVar./seedCntsThr(:, 2).^2));

                pvalSpecs=[lgnMu, repelem(lgnSigma(1), length(lgnMu),1)];

            else
                seedCntsThr(:, 2)=seedCntsThr(:, 2)*(options.gTrainNum+addPos);
                pvalSpecs=[nPos, nNeg];
            end
            pvalues=computePvalue(seedCntsThr, pvalSpecs, bernoulli);
            else
                pvalues=ones(size(appxPvalues));
            end

        else
            pvalues=appxPvalues;
        end

        if options.diffMotif
            pvaluesE=[pvalues;100];
            pvaluesd=pvaluesE(2:end)-pvaluesE(1:end-1);
            pvalInc=options.pvalIncRatio*pvaluesE(1:end-1);
            pvalInc(pvalInc>-1)=-1;

            idx=find(pvaluesd>pvalInc, 1);
            pvaluemin=pvalues(idx);
        else
            [pvaluemin, idx]=min(pvalues);
        end
        pvalue0=pvalue;
        pvalue=pvaluemin;
        thrOptimum0=thrOptimum;
        thrOptimum=PWMSTV(end-idx+1);
        selIndexes=(appxPWMS>=thrOptimum);
        sumSelIndexes=sum(selIndexes);
        appxSeeds=appxSeeds(1:sumSelIndexes, :);
        seedsPvalue0=appxSeedsLogPvalues;
        appxSeedsLogPvalues=appxPvalues(1:sumSelIndexes);
        appxSeedsLogPvalues=min(appxSeedsLogPvalues, 0);
        seedsToPWM0=seedsToPWM;
        PWMOut=PWM1;
        if (pvalue<pvalue0 && refIter<nRefIter)
            PWM10=PWM1;
            if sumSelIndexes>1
                zoopsCnt=zoopsCntAll(1:sumSelIndexes,1);
                totalCnt=totalCntAll(1:sumSelIndexes,1);
                totalCnt(totalCnt<1)=1;
                wgtPosCount=-zoopsCnt./totalCnt.*appxSeedsLogPvalues;
                seedsToPWM=appxSeeds;
                numSelIndex=size(seedsToPWM,1);

                seedCellTi=seedsToPWM(:)==cellCharsDouble;
                PWM1T=seedCellTi.*repmat(wgtPosCount(:), wMax,1);
                PWM1T=reshape(PWM1T,numSelIndex,wMax, numCells);
                PWM1=squeeze(sum(PWM1T, 1));
                PWM1=PWM1.';
                PWM1=PWM1+prior*PWM0;
                PWM1=PWM1./sum(PWM1);

                if (options.isPal)
                    PWM1=(PWM1+PWM1(:, end:-1:1))/2;
                end



                PWM1L=log2(PWM1);
                if options.mkvOrder>0
                    PWMS=PWM1L;
                else
                    PWMS=(PWM1L-PWM0L);
                end


            else
                seedsToPWM=appxSeeds;
                wgtPosCount=-appxSeedsLogPvalues;
                PWM1=wgtPosCount*(seedsToPWM(:)==cellCharsDouble);
                PWM1=PWM1.';
                PWM1=PWM1+prior*PWM0;
                PWM1=PWM1./sum(PWM1);

                if (options.isPal)
                    PWM1=(PWM1+PWM1(:, end:-1:1))/2;
                end

                PWM1L=log2(PWM1);
                if options.mkvOrder>0
                    PWMS=PWM1L;
                else
                    PWMS=(PWM1L-PWM0L);
                end
            end
            refIter0=refIter;
            refIter=refIter+1;
        else
            refIter0=refIter;

            refIter=nRefIter+1;

        end


    end

end

if refIter0==1
    [~, consensusSeed]=max(PWM1);
    seedEnriched.seed=consensusSeed;
    iterNumber=refIter;
    seedEnriched.pvalue=pvalue;
    seedEnriched.PWMS=PWM1;
    seedEnriched.thresholdOptimum=thrOptimum;
    seedEnriched.seedsToPWM=seedsToPWM;
    seedEnriched.PWMOut=PWM1;
else

    [~, consensusSeed]=max(PWM10);
    seedEnriched.seed=consensusSeed;
    iterNumber=refIter0;
    seedEnriched.pvalue=pvalue0;
    seedEnriched.PWMS=PWM10;
    seedEnriched.thresholdOptimum=thrOptimum0;
    seedEnriched.seedsToPWM=seedsToPWM0;
    seedEnriched.PWMOut=PWMOut;

end
seedEnriched.seedsPvalue=seedsPvalue0;
end
