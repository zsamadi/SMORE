function outMotif=scoreModelRobust(seedsEnriched, seqData, options)


options.isUHold=true;

if ~isempty(seedsEnriched.seeds)

    nPos=floor(options.numHNodes);
    nNeg=nPos;

    options.trNumShuffle=1;

    shConfig=options.shConfig;


    shConfig.fixedNodes=[];
    shConfig.numShuffle=1;
    addPos=options.addPos;

    seedsToPWM=seedsEnriched.seedsToPWM{1};


    scrSpecs.back=seqData.back;
    scrSpecs.mkvOrder=seqData.mkvOrder;
    scrSpecs.rvp=options.rvp;
    scrSpecs.seedsToPWM=seedsToPWM;
    PWMSCell=seedsEnriched.PWMSCell;
    PWMSOutCell=seedsEnriched.PWMSOutCell;
    thrOpt=seedsEnriched.thresholdOptimum;
    if nPos>0 && options.scIterMax>0

        if scrSpecs.mkvOrder>0
            PWMSi=log2(PWMSCell{1});
        else
            PWMSi=log2(PWMSCell{1})-log2(scrSpecs.back{1});
        end

        options.PWMS=PWMSi;
        options.isPos=true;
        cPHTypes=seqData.cPHTypes;

        unqMersP=getPSNmers(seqData, cPHTypes,options);

        options.isPos=false;

        [posCntN, pSeedsInPWM]=siteCountRobust(unqMersP,seedsToPWM);
        outMotif.testPvalueVec=zeros(options.scIterMax, 1);

        pvalReg=zeros(options.scIterMax,1);
        nSeedsInPWMReg=cell(options.scIterMax,1);
        pvalSpecReg=zeros(options.scIterMax,2);


        seqData.cNHTypes=getShuffleTypes(seqData.cNHTypes, shConfig);

        

        for iterSC=1:options.scIterMax
            negCntNM=zeros(options.gTrainNum+addPos, 1);
            nSeedsInPWM0C=cell(options.gTrainNum+addPos, 1);
            fprintf('%3d%%',round(iterSC/options.scIterMax*100));

            negCntNM(1)=posCntN;
            nSeedsInPWM0C{1}=pSeedsInPWM(:, end);



            for iterGT=1:options.gTrainNum
                cNHTypes=getShuffleTypes(cPHTypes, shConfig);
                unqMersN=getPSNmers(seqData,cNHTypes, options);


                [negCntN, nSeedsInPWM0]=siteCountRobust(unqMersN, seedsToPWM);

                negCntNM(iterGT+addPos)=negCntN;
                nSeedsInPWM0C{iterGT+addPos}=nSeedsInPWM0(:, end);

                
            end




                bernoulli0=options.bernoulli;
                pvalSpecs=[nPos, nNeg];


                pvalSpecReg(iterSC, :)=[bernoulli0, 0];
                negCntN=sum(negCntNM);



              if posCntN>0
                    PEnrich=computePvalue([posCntN, negCntN], pvalSpecs, bernoulli0);
                else
                    PEnrich=0;
               end



            pvalReg(iterSC)=PEnrich;
            nSeedsInPWM0=horzcat(nSeedsInPWM0C{:});

            uNMers=seedsToPWM;

            nWeights=sum(nSeedsInPWM0, 2);


            nSeedsInPWMReg{iterSC}=[uNMers, nWeights];
            fprintf('\b\b\b\b')

        end
        [pvalReg, istReg]=sort(pvalReg);
        pvalSpecReg=pvalSpecReg(istReg, :);
        nSeedsInPWMReg=nSeedsInPWMReg(istReg);
        if options.scIterMax>=10
            perceID=floor(0.95*options.scIterMax);
            perceID=max(perceID, 1);
        else
            perceID=options.scIterMax;
        end

        % outMotif.testPvalue=max(pvalReg(perceID), seedsEnriched.pvalues(1));
        outMotif.testPvalue=pvalReg(perceID);
        nSeedsInPWM=nSeedsInPWMReg{perceID};
        pvalSpecs=pvalSpecReg(perceID, :);

        outMotif.testPvalueVec=pvalReg;

        pnCounts=zeros(size(seedsToPWM, 1), 3);

        pnCounts(:, 1)=pSeedsInPWM(:, end);


        pnCounts(:, 2)=nSeedsInPWM(:, end);

        if options.isNormalTest
            indvPvalues=1-logncdf(pnCounts(:, 1),pvalSpecs(:, 1),pvalSpecs(:, 2));
        else
            indvPvalues=computePvalue(pnCounts(:, 1:2), [nPos, nNeg], options.bernoulli);
        end
        pnCounts(:, end)=indvPvalues;
        outMotif.pnsCPV=pnCounts;


    else
        PEnrich=seedsEnriched.pvalues(1);
        outMotif.testPvalue=PEnrich;
        outMotif.testPvalueVec=PEnrich;
        outMotif.pnsCPV=[0,0];
    end

    outMotif.PWMSE=log2(PWMSCell{1})-log2(scrSpecs.back{1});
    outMotif.cSeed=seedsEnriched.seeds(1, :);
    outMotif.scoreThr=thrOpt(1);
    outMotif.trainPvalue=seedsEnriched.pvalues(1);
    outMotif.PWMSOut=log2(PWMSOutCell{1})-log2(scrSpecs.back{1});
    outMotif.secSeeds=seedsEnriched.seeds(2:end, :);

    outMotif.seedsToPWM=seedsEnriched.seedsToPWM{1};
    outMotif.seedsPvalue=seedsEnriched.seedsPvalue{1};
else

    outMotif=[];
end