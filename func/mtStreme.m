function [outMotifCell,textOut,commandText, background, cPTypesOut]= mtStreme(seqData, options)

% main function for applying streme on the input data
% Input arguments:
%   - seqData:
%       - seqData.pSeq  : primary sequence
%       - seqData.nSeq  : control sequence
%       - seqData.pHSeq : primary holdoput sequence
%       - seqData.nHSeq : control holdoput sequence
%       - seqData.lens  : data lengths [size(posSeq, 1),size(negSeq, 1),size(pHoldSeq, 1),size(nHoldSeq, 1)];
%   - options
%       - options.NEVAL     : Number of motifs for initial evaluation
%       - options.NREF      : Number of motifs for refinement
%       - options.prior     : Dirichlet prior for initial evaluation
%       - options.REFPrior  : Dirichlet prior for refinement step
% Output arguments:
%   - outMotifCell: Cell of outMotifs
%       - outMotif.PWMSE        : Final refined PWM score matrix
%       - outMotif.consensusSeed: Final consensus Seed
%       - outMotif.scoreThr     : score threshold used for erasing
%       - outMotif.testPvalue   : significance obtained from hold-out
%       - outMotif.trainPvalue  : significance obtained from training
%       - textOut               : output text indicating the reason for finishing the algorithm
%       - commandText               : the command that generated the output results


arguments
    seqData.cPTypes  (:, 1) double
    seqData.cPHTypes  (:, 1) double
    seqData.cNTypes  (:, 1) double
    seqData.cNHTypes  (:, 1) double
    seqData.pSeq (:, 1) cell
    seqData.posWeight=[];
    seqData.negWeight=[];
    seqData.pHoldWeight=[];
    seqData.nHoldWeight=[];
    seqData.fixedNodes=[];
    options.alphabet='ABCDEFGHIJKLMNO';
    options.NEVAL=25;
    options.NREF=4;
    options.prior=0.01;
    options.REFPrior=0.1;
    options.ENRCHPrior=0.1;
    options.wMin=3;
    options.wMax=3;
    options.MinSeedWidth=3;
    options.nRefIter=20;
    options.threshold=0.05;
    options.patience=3;
    options.evalue=false;
    options.nmotifs=0;
    options.zeropad=true;
    options.sortmwidth=true;
    options.eps=1e-16;
    options.numPartitions=2;
    options.mkvOrder=0;
    options.rvp=false;
    options.maxPalED=5;
    options.minPalRatio=10;
    options.isPal=false;
    options.bernoulli=-1;
    options.isFastZNIC=false;
    options.isErzNegative=true;
    options.isErzHOut=true;

    options.scIterMax=2;
    options.gTrainNum=50;
    options.isEraseNodes=true;
    options.fixedTypes=1;
    options.diffMotif=false;
    options.fileID=1;
    options.shConfig=[];
    options.shuffleMode='shuffle';
    options.isKernel=false;
    options.isOLess=false;
    options.isEnrich=true;
    options.isErzFixNodes=false;
    options.isUBack=false;
    options.trNumShuffle=1;
    options.pvalIncRatio=0.01;
    options.indSeedMode=false;
    options.isOKSingleNC=true;
    options.isUHold=false;
    options.isNormalTest=true;

end

if strcmpi(options.shuffleMode, 'kernel')
    options.isKernel=true;
end

betaP=1;
alphaP=10;

% logFileID=options.fileID;

options.MinSeedWidth=min(options.wMin, options.MinSeedWidth);
options.addPos=1;

% options.MinSeedWidth=options.wMin;

%
%     pSeqName=inputname(2);
%
%     nSeqName=inputname(4);
%
%      pHSeqName=inputname(6);
%
%     nHSeqName=inputname(8);
%
%     PWM0Name=inputname(10);
%     alphabetName=inputname(18);


% commandText= sprintf("mStreme(pSeq=%s ,nSeq=%s,pHSeq=%s, nHSeq=%s, PWM0=%s, wMin=%d, wMax=%d,threshold=%2.2f, alphabet=%s)",...
%     pSeqName, nSeqName,pHSeqName, nHSeqName,PWM0Name, options.wMin, options.wMax, options.threshold,alphabetName);


% [posSeq, negSeq, pHoldSeq, nHoldSeq]=generateSeqs(seqData.pSeq, hofract);



commandText= sprintf("mStreme(pSeq={input data}  ,nSeq=={shuffled data},pHSeq={positive holdout seq}, nHSeq={control holdout seq},...\n PWM0={bg frequency}, wMin=%d, wMax=%d,threshold=%10.10f, alphabet={input alphabet})",...
    options.wMin, options.wMax, options.threshold);



background=getMarkovFromSequence(seqData.pSeq,seqData.cPTypes,options.alphabet, options.mkvOrder);

if options.isUBack
    tmp=background{1};
    tmp=ones(size(tmp));
    tmp=tmp/sum(tmp);
    background{1}=tmp;
end

options.numCells=length(options.alphabet);
options.numNodes=length(seqData.cPTypes);

if options.rvp
    numPSeqs=length(seqData.pSeq)/2;
    numNSeqs=numPSeqs;
    numPHSeqs=numPSeqs;
    numNHSeqs=numPSeqs;
else
    numPSeqs=length(seqData.pSeq);
    numNSeqs=numPSeqs;
    numPHSeqs=numPSeqs;
    numNHSeqs=numPSeqs;
end

seqData.lens=[numPSeqs,numNSeqs,numPHSeqs,numNHSeqs];
options.numHNodes=numPHSeqs;

seqData.mkvOrder=options.mkvOrder;
seqData.back=background;
seqData.rvp=options.rvp;


wMin=options.wMin;
wMax=options.wMax;

isEraseOccured=true;
iErase=1;
outMotifCell=cell(100,wMax-wMin+1);
% cSeeds=zeros(100, wMax);
% trainPvalues=zeros(100, wMax-wMin+1);
testPvalues=zeros(100, wMax-wMin+1);


validThreshold= log(options.threshold);
numPateince=0;


% coptions.numCells=options.numCells;
% coptions.wMax=options.wMax;
% coptions.MinSeedWidth=options.MinSeedWidth;
% coptions.numSeqs=numPSeqs;
% coptions.trNumShuffle=options.trNumShuffle;
% coptions.rvp=options.rvp;

% coptions.numNodes=length(seqData.cPTypes);
% coptions.isFastZNIC=options.isFastZNIC;

% coptions.isProbCount=false;

% coptions.rvp=false;

%
% coptions.modCount=false;
% [unqMersCell, wMerZoops]=countSeeds(seqData, coptions);




% load('countSeeds.mat','unqMersCell', 'wMerZoops')



seqDataSpecs.back=background;
seqDataSpecs.lens=seqData.lens;
seqDataSpecs.mkvOrder=seqData.mkvOrder;


% eznConfigs.rvp=options.rvp;

% eznConfigs.back=background;

% unqMersFlagsCellIn=cell(wMax-options.MinSeedWidth+1,1);
% coptions.isPos=true;
% coptions.isPN=false;
% unqMersP=countSeedIs(seqData,cell(wMax-options.MinSeedWidth+1,1), coptions);

% coptions.isPos=false;
% unqMersN=countSeedIs(seqData,unqMersP, coptions);

% for iLt=1:length(unqMersP)
%     seedsWmaxNi.seeds=unqMersP{iLt}.seeds;
%     seedsN=unqMersN{iLt}.seeds;
%     [sFlags, iSNg]=ismember(unqMersP{iLt}.seeds,seedsN, 'rows');
%     nCounts=zeros(size(sFlags));
%     nCounts(sFlags)=unqMersN{iLt}.counts(iSNg, 1);
%     seedsWmaxNi.counts=[unqMersP{iLt}.counts(:, 1), nCounts];



% coptions.isPN=true;
% coptions.isOLess=options.isOLess;

% coptions.trNumShuffle=options.trNumShuffle;


% coptions.cPTypesInit=seqData.cPTypes;
% coptions.cNTypesInit=seqData.cNTypes;

options.cPTypesInit=seqData.cPTypes;
options.cNTypesInit=seqData.cNTypes;

options.cPTypesInit=seqData.cPTypes;
options.cNTypesInit=seqData.cNTypes;
% coptions.isOKSingleNC=options.isOKSingleNC;
isTrainGraph=false;


if isTrainGraph

    [unqMersP, unqMersNC, weightsM]=trainGraph(seqData, options);

    contactTypeMatrix=zeros(length(background{1}));

    for iSeed=1:length(unqMersP.seeds)
        contactTypeMatrix(unqMersP.seeds(iSeed, 1), unqMersP.seeds(iSeed, 2))=unqMersP.weights(iSeed, 1);
    end


    contactTypeMatrix=contactTypeMatrix./sum(contactTypeMatrix, 2)*100;

    contactTypeMatrix(contactTypeMatrix<1)=0;
    contactTypeMatrix=contactTypeMatrix./sum(contactTypeMatrix, 2)*100;


    % figure
    % heatmap(contactTypeMatrix)
    % clim([0 50])
    % colormap(flipud(gray))
    %
    %
    % load('znicCNTMap.mat')
    %
    % contactTypeMatrixDiff=contactTypeMatrix-znicCNTMap;
    %
    %
    % figure
    % heatmap(contactTypeMatrixDiff)
    % colormap("jet")
    % title('All minus ZNIC')


    weightsMean=mean(weightsM, 2);
    weightsSTD=std(weightsM,0, 2);
end





options.isNotNMer=true;

cNTypesMat=zeros(length(seqData.cPTypes),options.gTrainNum);

for iShuffle=1:options.gTrainNum

    cNTypesi=getShuffleTypes(seqData.cPTypes, options.shConfig);
    cNTypesMat(:,iShuffle)=cNTypesi;
end

options.cNTypesMat=cNTypesMat;

fprintf('\n')
while(isEraseOccured && numPateince<options.patience)

    for iW=1:wMax-wMin+1

        fprintf('Counting Seeds         (Motif #%d)... ', iErase);
        ticCS=tic;


        % seedsWmax=countSeedIs(seqData,unqMersP, coptions);

        [seedsWmax, ~, ~, seqData]=countPNSeeds(seqData, options);
        options.isNotNMer=false;



        % [seedsWmax, ~]=countSeeds(seqData,unqMersFlagsCellIn, coptions);
        %         unqMersFlagsCellIn=cell(wMax-options.MinSeedWidth+1,1);
        elapsedTime=toc(ticCS);
        fprintf('Elapsed Time %3.3f seconds \n' , elapsedTime);




        % Finding first NEVAL significant seeds
        %
        if iErase==1 || options.isErzFixNodes

            totalScale=sum(seedsWmax{end}.weights, 1);
            % totalScale(2)=totalScale(2)*(options.gTrainNum+1);
            options.posTotalC=totalScale;
            % backFreq=background{1};

            % backFreq=backFreq*length(seqData.cPTypes);

            % posCounts=seedsWmax{end}.weights(:, 1);

            % validSeeds=posCounts>0;

            % posCounts=posCounts(validSeeds);



            % allSeeds=seedsWmax{end}.seeds;
            % allSeeds=allSeeds(validSeeds, :);

            % allSeedsFreq=backFreq(allSeeds);
            % allSeedsFreq=prod(allSeedsFreq, 2);

            % allSeedsFreq=allSeedsFreq./min(allSeedsFreq);

            % betaN=betaP+sum(allSeedsFreq);           

            % ynp1=seedsWmax{end}.weights(:, 1);
            % 
            % alphaN=alphaP+totalScale(2);
            % PBayes=betaN./(betaN+allSeedsFreq);


            % ynp1FromAphanp1=nchoosekLog(alphaNp1-1, ynp1);

            % pYnp1=ynp1FromAphanp1+alphaN*betaNLog-alphaNp1*betaNp1Log;
            % 
            % pYnp1 = nbincdf(ynp1,alphaN,PBayes,'upper');
            % pYnp1=log(pYnp1);






            % posCountsPerPopl=posCounts./allSeedsFreq;

            % alphaOBeta=mean(posCountsPerPopl);

            % oneOPopuMean=mean(1./allSeedsFreq);
            % varYOPopu=var(posCountsPerPopl);

            % alphaOBeta2=varYOPopu-oneOPopuMean*alphaOBeta;
            % betaE=alphaOBeta/alphaOBeta2;

            % alphaE=alphaOBeta*betaE;

            % thetaMeanE=(alphaE+posCounts)./(betaE+allSeedsFreq);

            % [thetaMeanES, iESt]=sort(thetaMeanE, 'descend');

            % allSeedSorted=allSeeds(iESt, :);

            % posCountSort=posCounts(iESt, :);

            % pd=fitdist(posCountsPerPopl, 'Gamma');




            check=1;




           








            


            % countsMat=unqMersCell{end}.counts;
            if options.isNormalTest
                options.bernoulli=2;
            else


                % options.bernoulli=[alphaE, betaE];
                options.bernoulli=totalScale(1)/sum(totalScale);

                % options.bernoulli=(1+options.gTrainNum)/(2+options.gTrainNum);
            end
        end

        % countsMat=unqMersCell{end}.counts;
        % totalScale=sum(seedsWmax{end}.counts);
        %
        % eoptions.bernoulli=totalScale(1)/sum(totalScale);
        % options.bernoulli=eoptions.bernoulli;


        fprintf('Evaluate Initial Seeds (Motif #%d)... ', iErase);


        ticEV=tic;


        if ~isempty(seedsWmax{end}.seeds)

            % [seedsNEVALB, seedsWmaxB, pvaluesRawOutB]=bayesInitialSeeds(seedsWmax, options);


            [seedsNEVAL, seedsWmax]=evaluateInitialSeeds(seedsWmax, options);

            % ind2460=find(all(sort(seedsNEVAL.seeds, 2)==[24, 60], 2));
            % if ~isempty(ind2460)
            % fprintf('*****************2460 in ierase %d is %d*********************\n', iErase, ind2460);
            % end

            elapsedTime=toc(ticEV);
            fprintf('Elapsed Time %3.3f seconds \n' , elapsedTime);


            % refining first NEVAL seeds to find first NREF seeds
            %         tic
            fprintf('Refine Initial Seeds   (Motif #%d)... ', iErase);

            ticRI=tic;
            seedsEnriched=refineInitialSeeds(seedsNEVAL,seedsWmax,seqData, seqDataSpecs, options);
            elapsedTime=toc(ticRI);
            fprintf('Elapsed Time %3.3f seconds \n' , elapsedTime);
            %         toc



            % Further refine first NREF seeds
            % uses following data from seqData
            %  PWM0=seqData.PWM0;
            %  nPos=seqData.lens(1);
            %  nNeg=seqData.lens(2);
            fprintf('Nested Subset Enrich   (Motif #%d)... ', iErase);
            %        if iErase==9
            %            check=1;
            %        end

            ticNS=tic;
            seedsEnriched=nestedSubsetEnrichment(seedsEnriched,seedsWmax,seqData,seqDataSpecs,options);
            elapsedTime=toc(ticNS);
            fprintf('Elapsed Time %3.3f seconds \n' , elapsedTime);


            % Enrichment test on first NREF motifs
            fprintf('Hold-out Scoring       (Motif #%d)... ', iErase);

            ticHS=tic;
            outMotif=scoreModelRobust(seedsEnriched, seqData, options);
            %         outMotif=scoreModelPssm(seedsEnriched, seqData, options);
            elapsedTime=toc(ticHS);
            fprintf('Elapsed Time %3.3f seconds \n' , elapsedTime);


            fprintf('Erase Found Motif      (Motif #%d)... ', iErase);
            ticEM=tic;
            [seqData,~, isEraseOccured]=eraseMotif(outMotif,  seqData, options);
            elapsedTime=toc(ticEM);
            fprintf('Elapsed Time %3.3f seconds \n' , elapsedTime);

            %         eznConfigs.cTypes=cTypesEx;



            %         nMersCellO=eraseMotifMers(outMotif,  nMersCell, eznConfigs);

            % outMotif.pNodes=sitesErased.pNodes;
            % outMotif.nNodes=sitesErased.nNodes;
            % outMotif.pHNodes=sitesErased.pHNodes;
            % outMotif.nHNodes=sitesErased.nHNodes;





            % outMotif.nsites=sum(sitesErased.pSeq)+sum(sitesErased.pHSeq);

            PWMSE=outMotif.PWMSOut;
            if options.mkvOrder>0
                PWM1=2.^(PWMSE);
            else
                PWM1=2.^(PWMSE).*background{1};
            end

            % PWM1=round(PWM1, 6);

            outMotif = rmfield(outMotif,'PWMSE');

            outMotif.PWM=PWM1;
            outMotif.iErase=iErase;


            outMotifCell{iErase, iW}=outMotif;
            % cSeeds(iErase, :)=outMotif.cSeed;
            % trainPvalues(iErase, :)=outMotif.trainPvalue;
            testPvalues(iErase, :)=outMotif.testPvalue;

            % Check if erasure has been occured

            % isEraseOccured=any(sitesErased.pSeq)||any(sitesErased.nSeq);
            % isEraseOccured=true;

            %         if ~isEraseOccured
            %             check=1;
            %         end

            %         if isEraseOccured
            %         coptions.modCount=true;
            %
            %             unqMersCellN=modifyCountSeeds(seqDataCopy, seqData,sitesErased,unqMersCell, coptions);
            %         end


            if (options.evalue)
                validThreshold= log(options.threshold)-log(iErase);
            end

            if options.nmotifs==0

                if outMotif.testPvalue>validThreshold
                    numPateince=numPateince+1;
                else
                    numPateince=0;
                end
            else
                if (iErase==options.nmotifs)
                    numPateince=options.patience +1;
                end
            end

        else
            isEraseOccured=false;
        end


    end

    % if options.isEraseNodes
    %     ppHNode=[outMotif.pNodes;outMotif.pHNodes];
    %     ppHNode=ppHNode(:);
    %     ppHNode=unique(ppHNode);
    %     seqData.cPTypes(ppHNode)=0;
    %     % seqData.cPHTypes(ppHNode)=0;
    %
    %     nnHNode=[outMotif.nNodes;outMotif.nHNodes];
    %     nnHNode=nnHNode(:);
    %     nnHNode=unique(nnHNode);
    %     seqData.cNTypes(nnHNode)=0;
    %     % seqData.cNHTypes(nnHNode)=0;
    %
    % end





    iErase=iErase+1;






end


cPTypesOut=seqData.cPTypes;


if ~isEraseOccured
    textOut='Motif Analysis Finished Becasue No Site Erased\n';
    iErase=iErase-1;
else

    if options.nmotifs>0
        textOut=sprintf('Motif Analysis Finished Becasue %d Motifs Were Found\n', options.nmotifs);

    else

        textOut=sprintf('Motif Analysis Finished Becasue %d Consecutive Nonsignificant Motifs Were Found\n', options.patience);

    end



end
fprintf(textOut);
mumFoundMotifs=iErase-1;
% cSeeds=cSeeds(1:mumFoundMotifs, :);
% trainPvalues=trainPvalues(1:mumFoundMotifs, :);
testPvalues=testPvalues(1:mumFoundMotifs, :);

outMotifCell=outMotifCell(1:iErase-1, :);
% seedsNEVALO=seedsNEVALO(1:iErase-1, :);
% validMotifs=testPvalues<=validThreshold;
%
% validMotifs=testPvalues<=0;
%
% testPvalues=testPvalues(validMotifs);
% outMotifCell=outMotifCell(validMotifs);
[~, idx]=sort(testPvalues);

% idx=(1:length(testPvalues));

outMotifCell=outMotifCell(idx);
% trainPvalues=trainPvalues(idx);
% consensusSeeds=consensusSeeds(idx, :);

% mumFoundMotifs=length(testPvalues);

% for im=1:mumFoundMotifs
%
%     outMotifCell{im}.testPvalue=log(mumFoundMotifs)+outMotifCell{im}.testPvalue;
% end
% check=1;

end



