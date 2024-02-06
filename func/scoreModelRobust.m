function outMotif=scoreModelRobust(seedsEnriched, seqData, options)

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

options.isUHold=true;

if ~isempty(seedsEnriched.seeds)

    nPos=floor(options.numHNodes);
    nNeg=nPos;

    % G=options.G;
    options.trNumShuffle=1;

    shConfig=options.shConfig;

    % cSections=ones(size(seqData.cPHTypes, 1), 1);
    % xyCoordinates=G.Nodes.Coordinates;
    % shConfig.cSections=cSections;

    shConfig.fixedNodes=[];
    shConfig.numShuffle=1;
    addPos=options.addPos;

    % if strcmpi(shConfig.shuffleMode, 'kernelPath')
    %     addPos=0;
    % else
    %     addPos=1;
    % end
    seedsToPWM=seedsEnriched.seedsToPWM{1};


    scrSpecs.back=seqData.back;
    scrSpecs.mkvOrder=seqData.mkvOrder;
    scrSpecs.rvp=options.rvp;
    scrSpecs.seedsToPWM=seedsToPWM;
    PWMSCell=seedsEnriched.PWMSCell;
    PWMSOutCell=seedsEnriched.PWMSOutCell;
    thrOpt=seedsEnriched.thresholdOptimum;
    if nPos>0 && options.scIterMax>1

        if scrSpecs.mkvOrder>0
            PWMSi=log2(PWMSCell{1});
        else
            PWMSi=log2(PWMSCell{1})-log2(scrSpecs.back{1});
        end

        options.PWMS=PWMSi;
        options.isPos=true;
        cPHTypes=seqData.cPHTypes;

        unqMersP=getPSNmers(seqData, cPHTypes,options);






        % options.isPos=true;
        % options.isPN=false;

        % unqMersP=countHPSeeds(seqData, options);
        options.isPos=false;
        % totalScale=zeros(1,2);
        % totalScale(1)=sum(unqMersP.counts(:, 1));

        totalScale=options.posTotalC;
        [posCntN, pSeedsInPWM]=siteCountRobust(unqMersP,seedsToPWM);
        outMotif.testPvalueVec=zeros(options.scIterMax, 1);

        pvalReg=zeros(options.scIterMax,1);
        nSeedsInPWMReg=cell(options.scIterMax,1);
        pvalSpecReg=zeros(options.scIterMax,2);


        seqData.cNHTypes=getShuffleTypes(seqData.cNHTypes, shConfig);

        
        % fprintf('score iter    ... ');


        for iterSC=1:options.scIterMax
            negCntNM=zeros(options.gTrainNum+addPos, 1);
            nSeedsInPWM0C=cell(options.gTrainNum+addPos, 1);
            fprintf('%3d%%',round(iterSC/options.scIterMax*100));

            negCntNM(1)=posCntN;
            nSeedsInPWM0C{1}=pSeedsInPWM(:, end);



            for iterGT=1:options.gTrainNum
                cNHTypes=getShuffleTypes(cPHTypes, shConfig);
                unqMersN=getPSNmers(seqData,cNHTypes, options);

                % unqMersN=countHPSeeds(seqData, options);


                [negCntN, nSeedsInPWM0]=siteCountRobust(unqMersN, seedsToPWM);

                negCntNM(iterGT+addPos)=negCntN;
                nSeedsInPWM0C{iterGT+addPos}=nSeedsInPWM0(:, end);

                
            end
             % totalScale(2)=totalScale(1);



            if options.isNormalTest
                negCntNMean2=(mean(negCntNM)).^2;
                negCntNVar=var(negCntNM);
                vm2=negCntNVar+negCntNMean2;


                lgnMu=log(negCntNMean2./sqrt(vm2));
                lgnSigma=sqrt(log(vm2./negCntNMean2));

                % seedCntsThrVar=seedCntsThrVar-(seedCntsThr(:, 2)).^2;

                % PEnrich= 1-logncdf(posCntN,lgnMu,lgnSigma);
                pvalSpecs=[lgnMu,lgnSigma];

                pvalSpecReg(iterSC, :)=[lgnMu,lgnSigma];
                bernoulli0=2;
                negCntN=mean(negCntNM);
            else

                bernoulli0=options.bernoulli;
                pvalSpecs=[nPos, nNeg];


                pvalSpecReg(iterSC, :)=[bernoulli0, 0];
                negCntN=sum(negCntNM);

            end

              if posCntN>0
                    PEnrich=computePvalue([posCntN, negCntN], pvalSpecs, bernoulli0);
                else
                    PEnrich=0;
               end



            pvalReg(iterSC)=PEnrich;
            nSeedsInPWM0=horzcat(nSeedsInPWM0C{:});

            uNMers=seedsToPWM;

            nWeights=sum(nSeedsInPWM0, 2);
            % nWeights=nWeights/(options.gTrainNum+addPos);


            nSeedsInPWMReg{iterSC}=[uNMers, nWeights];
            fprintf('\b\b\b\b')

        end
        % fprintf('\n');

        [pvalReg, istReg]=sort(pvalReg);
        pvalSpecReg=pvalSpecReg(istReg, :);
        nSeedsInPWMReg=nSeedsInPWMReg(istReg);

        perceID=floor(0.95*options.scIterMax);

        outMotif.testPvalue=pvalReg(perceID);
        nSeedsInPWM=nSeedsInPWMReg{perceID};
        pvalSpecs=pvalSpecReg(perceID, :);

        outMotif.testPvalueVec=pvalReg;

        pnCounts=zeros(size(seedsToPWM, 1), 3);

        % seedsToPWMSt=scrSpecs.seedsToPWM;
        % pSeedsInPWMSt=pSeedsInPWM;
        % pILa=ismember(seedsToPWMSt, pSeedsInPWMSt(:, 1:end-1), 'rows');
        pnCounts(:, 1)=pSeedsInPWM(:, end);

        % nSeedsInPWMSt=nSeedsInPWM;

        pnCounts(:, 2)=nSeedsInPWM(:, end);

        if options.isNormalTest
            indvPvalues=1-logncdf(pnCounts(:, 1),pvalSpecs(:, 1),pvalSpecs(:, 2));
        else
            indvPvalues=computePvalue(pnCounts(:, 1:2), [nPos, nNeg], options.bernoulli);
        end
        pnCounts(:, end)=indvPvalues;
        outMotif.pnsCPV=pnCounts;

        % if any(all(pnCounts==0, 2))
        %     check=1;
        % end


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