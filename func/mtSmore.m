function [outMotifCell,textOut,commandText, background, cPTypesOut]= mtSmore(seqData, options)



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
    options.REFPrior=0.01;
    options.ENRCHPrior=0.01;
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
    options.fixedTypes=0;
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
    options.isNormalTest=false;
    options.iSyn=false;

end

if strcmpi(options.shuffleMode, 'kernel')
    options.isKernel=true;
end

options.MinSeedWidth=min(options.wMin, options.MinSeedWidth);
options.addPos=1;

commandText= sprintf("smore(pSeq={input data}  ,nSeq=={shuffled data},pHSeq={positive holdout seq}, nHSeq={control holdout seq},...\n PWM0={bg frequency}, wMin=%d, wMax=%d,threshold=%10.10f, alphabet={input alphabet})",...
    options.wMin, options.wMax, options.threshold);



background=getMarkovFromSequence(seqData.pSeq,seqData.cPTypes,options.alphabet, options.mkvOrder);

if options.isUBack

    tmp=background{1};
    tmp=ones(size(tmp));
    tmp=tmp/sum(tmp);
    background{1}=tmp;
elseif options.iSyn
    tmp=1./((1:length(background{1})).^(1/4));
    tmp=tmp/sum(tmp);
    tmp=reshape(tmp, size(background{1}));
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
testPvalues=zeros(100, wMax-wMin+1);
validThreshold= log(options.threshold);
numPateince=0;
seqDataSpecs.back=background;
seqDataSpecs.lens=seqData.lens;
seqDataSpecs.mkvOrder=seqData.mkvOrder;
options.cPTypesInit=seqData.cPTypes;
options.cNTypesInit=seqData.cNTypes;
options.cPTypesInit=seqData.cPTypes;
options.cNTypesInit=seqData.cNTypes;
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

        [seedsWmax, ~, ~, seqData]=countPNSeeds(seqData, options);
        options.isNotNMer=false;

        elapsedTime=toc(ticCS);
        fprintf('Elapsed Time %3.3f seconds \n' , elapsedTime);

        if iErase==1 || options.isErzFixNodes

            totalScale=sum(seedsWmax{end}.weights, 1);
            options.posTotalC=totalScale;


            if options.isNormalTest
                options.bernoulli=2;
            else

                options.bernoulli=totalScale(1)/sum(totalScale);
            end
        end



        fprintf('Evaluate Initial Seeds (Motif #%d)... ', iErase);


        ticEV=tic;


        if ~isempty(seedsWmax{end}.seeds)

            [seedsNEVAL, seedsWmax]=evaluateInitialSeeds(seedsWmax, options);

            elapsedTime=toc(ticEV);
            fprintf('Elapsed Time %3.3f seconds \n' , elapsedTime);
            fprintf('Refine Initial Seeds   (Motif #%d)... ', iErase);

            ticRI=tic;
            seedsEnriched=refineInitialSeeds(seedsNEVAL,seedsWmax,seqData, seqDataSpecs, options);
            elapsedTime=toc(ticRI);
            fprintf('Elapsed Time %3.3f seconds \n' , elapsedTime);
            fprintf('Nested Subset Enrich   (Motif #%d)... ', iErase);

            ticNS=tic;
            seedsEnriched=nestedSubsetEnrichment(seedsEnriched,seedsWmax,seqData,seqDataSpecs,options);
            elapsedTime=toc(ticNS);
            fprintf('Elapsed Time %3.3f seconds \n' , elapsedTime);

            fprintf('Hold-out Scoring       (Motif #%d)... ', iErase);

            ticHS=tic;
            outMotif=scoreModelRobust(seedsEnriched, seqData, options);
            elapsedTime=toc(ticHS);
            fprintf('Elapsed Time %3.3f seconds \n' , elapsedTime);


            fprintf('Erase Found Motif      (Motif #%d)... ', iErase);
            ticEM=tic;
            [seqData,~, isEraseOccured]=eraseMotif(outMotif,  seqData, options);
            elapsedTime=toc(ticEM);
            fprintf('Elapsed Time %3.3f seconds \n' , elapsedTime);
            PWMSE=outMotif.PWMSOut;
            if options.mkvOrder>0
                PWM1=2.^(PWMSE);
            else
                PWM1=2.^(PWMSE).*background{1};
            end

            outMotif = rmfield(outMotif,'PWMSE');

            outMotif.PWM=PWM1;
            outMotif.iErase=iErase;


            outMotifCell{iErase, iW}=outMotif;
            testPvalues(iErase, :)=outMotif.testPvalue;

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

testPvalues=testPvalues(1:mumFoundMotifs, :);

outMotifCell=outMotifCell(1:iErase-1, :);

[~, idx]=sort(testPvalues);

outMotifCell=outMotifCell(idx);


end



