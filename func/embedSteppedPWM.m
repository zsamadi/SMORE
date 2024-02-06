function [vertxType, embMotifC, PWMBgOut, patternCandOut]=embedSteppedPWM(pathList,numVertex, options)
config.mkvOrder=0;
node0=min(min(pathList));
pathList=pathList-node0+1;
numCells=options.nTypes;
[numPaths, W]=size(pathList);
numCandNodes=options.embedNumber;
patternNum=ceil(numCandNodes/W);


% cretae a matrix of all possible candiates for pwm embedding


allPatternCandTypesC=cell(W, 1);

randIDX=repmat((1:numCells).', 1, W);

for iff=1:W
    tmp=randIDX(:,iff);
    allPatternCandTypesC{iff}=(repmat(repelem(tmp, numCells^(W-iff)),numCells^(iff-1), 1)).' ;
end
allPatternCandTypes=vertcat(allPatternCandTypesC{:});
allPatternCandTypes=allPatternCandTypes.';
% allPatternCandTypes=allPatternCandTypes(allPatternCandTypes(:, 1)<=allPatternCandTypes(:, end), :);
allPatternCandTypes=unique(allPatternCandTypes, 'rows');

% first embed the  consensus seed

mtotifPattern=[5, 7, 8, 9, 2, 11, 12, 3, 1,6,  4, 10, 13, 14, 15];
patternRef=mtotifPattern(1:W);
PWMEpRef=zeros(W, numCells);
for iip=1:W 
    PWMEpRef(iip,patternRef(iip))=1;
end

% iteratively embed more deep pwms
maxDepth=options.maxDepth;
patternCandTypesC=cell(maxDepth, 1);
embMotifC=cell(maxDepth, 1);
patternCandTypesC{1}=repmat(patternRef, patternNum, 1);

patternCandTypesC1=cell(maxDepth, 1);
patternCandTypesC1{1}=patternRef;


prior=options.prior;
priorBg=prior/(1+prior)*options.PWM0;
PWMRef=PWMEpRef/(1+prior)+priorBg;    
PWMS=(log2(PWMRef)-log2(options.PWM0)).';


embMotif.PWM=PWMRef;
embMotif.cSeed=patternRef;
embMotif.scThr=8/10*scoreWords(patternRef, PWMS, config);
embMotif.bSeeds=patternRef;


embMotifC{1}=embMotif;
patternCandTypes0=patternRef;

minType=min(min(allPatternCandTypes));

allPatternCandTypesIndex=sum((allPatternCandTypes-minType).*(numCells.^(W-1:-1:0))+1, 2);
typeCnt=zeros(length(allPatternCandTypesIndex),1);
patternRefIndex=sum((patternRef-minType).*(numCells.^(W-1:-1:0))+1, 2);

typeCnt(allPatternCandTypesIndex==patternRefIndex)=patternNum;
patternNumi=patternNum;
for idepth=2:maxDepth

%     patternNumi=ceil((numCells-1)*0.3*patternNumi);

    allPatternCandTypesL=allPatternCandTypes;
    patternCandScores=scoreWords(allPatternCandTypesL, PWMS, config);    
    
    [patternCandScores, iSCR]=sort(patternCandScores, 'descend');
    patternCandScoresU=unique(patternCandScores);
    scrThr=patternCandScoresU(end-1);

    allPatternCandTypesL=allPatternCandTypesL(iSCR, :);

    
    allPatternCandTypesL=allPatternCandTypesL(patternCandScores>=scrThr, :);



    allPatternCandTypesL=setdiff(allPatternCandTypesL, patternCandTypes0, 'rows');
    allPatternCandTypesL=setdiff(allPatternCandTypesL, patternCandTypes0(:, end:-1:1), 'rows');
    allPatternCandTypesL=allPatternCandTypesL(8, :);
    patternCandTypesC1{idepth}=allPatternCandTypesL;

    numCandTypes=size(allPatternCandTypesL, 1);

%     patternNumi=ceil(patternNumi/numCandTypes);
    allPatternCandTypesL=repelem(allPatternCandTypesL, ceil(patternNumi/numCandTypes),1);




    % patternCandTypes(1, :)=[];
    
    
    
%     numCandTypes=size(allPatternCandTypesL, 1);
%     allPatternCandTypesL=allPatternCandTypesL(randperm(numCandTypes), :);
    allPatternCandTypesL=allPatternCandTypesL(1:patternNumi, :);
% 
% 
%     if numCandTypes==0
%         allPatternCandTypesL=allPatternCandTypesLPWM;
%     end
%     numCandTypes=size(allPatternCandTypesL, 1);


    
%     allPatternCandTypesL=allPatternCandTypesL(1:patternNumi, :);

    [allPatternCandTypesLU, ~, jU]=unique(allPatternCandTypesL, "rows");
    cntL=accumarray(jU, 1);
    allPatternCandTypesLUIndex=sum((allPatternCandTypesLU-minType).*(numCells.^(W-1:-1:0))+1, 2);

    cntIndex=ismember(allPatternCandTypesIndex, allPatternCandTypesLUIndex);

    typeCnt(cntIndex)=typeCnt(cntIndex)+cntL;

    patternCandTypesC{idepth}=allPatternCandTypesL;
    patternCandTypes0=vertcat(patternCandTypesC{:});
    patternCandTypes0=unique(patternCandTypes0, 'rows');



    allPatternCandTypesLPWMIndex=sum((patternCandTypes0-minType).*(numCells.^(W-1:-1:0))+1, 2);
    pwmIndex=ismember(allPatternCandTypesIndex, allPatternCandTypesLPWMIndex);

    allPatternCandTypesLN=repelem(patternCandTypes0, typeCnt(pwmIndex),1);
  
    embedNumberi=size(allPatternCandTypesLN, 1);
    
    patternCandTypes1=vertcat(patternCandTypesC1{:});
    
    PWM1T=patternCandTypes1(:)==(1:numCells);
    PWM1T=reshape(PWM1T,size(patternCandTypes1, 1),W, numCells);
    PWMPtOut=squeeze(sum(PWM1T, 1));
    %         PWMPtOut=PWMPtOut+options.prior*PWMBg0;
    PWMPtOut=PWMPtOut./sum(PWMPtOut,2);
    PWMPtOut=PWMPtOut/(1+prior)+priorBg;
    PWMS=(log2(PWMPtOut)-log2(options.PWM0)).';
    
    embMotif.PWM=PWMPtOut;
    embMotif.cSeed=patternRef;
    embMotif.scThr=scrThr;
    embMotif.bSeeds=patternCandTypes0;
    embMotifC{idepth}=embMotif;


end


patternCandTypes=vertcat(patternCandTypesC{:});
patternNum=size(patternCandTypes, 1);





shuffleIndex=randperm(numPaths);

pathListShuffle=pathList(shuffleIndex, :);
shuffleIndex=randperm(numPaths);
pathListShuffle=pathListShuffle(shuffleIndex, :);



if options.isBg
    patternNCand=pathListShuffle;
else



    patternCandIndex=zeros(patternNum, 1);
    
    patternCand=pathListShuffle;
    patternCandF=zeros(patternNum, W);
    patternCandF(1, :)=patternCand(1, :);
    iPCE=2;
    iPC=2;
    L0=W;
    patternCandIndex(1)=1;


    while(iPCE<=patternNum)
        
    
        nodesAll=[patternCandF(1:iPCE-1, :);patternCand(iPC, :)];
        nodesAll=unique(nodesAll(:));
        if length(nodesAll)==(L0+W)
            patternCandF(iPCE, :)=patternCand(iPC, :);
            patternCandIndex(iPCE)=iPC;
            iPCE=iPCE+1;
            L0=L0+W;
        end
        iPC=iPC+1;

    end

    patternCand=patternCandF;










    
    patternNCandIndex=(1:numPaths);
    patternNCandIndex(patternCandIndex)=[];
    patternNCand=pathListShuffle(patternNCandIndex, :);
end

vertxType=zeros(1, numVertex);

numAlph=options.nTypes;





PWMBg=options.PWM0;

errorNorm=1;

iterMax=100;
iter=1;
alphab=0.05;

if options.isBg
    errorNormVec=zeros(iterMax, 1);

    PWMBg0=PWMBg;
    
    while(errorNorm>0.001 && iter<iterMax)
        vertxType=embedPWM(patternNCand, PWMBg, vertxType, []);
        
        zeroType=(vertxType==0);
        numZeros=sum(zeroType);
    
        if numZeros>0
    
            vertxType(zeroType)=randi(numAlph, 1, numZeros);
        end

        patternNCandType=vertxType(patternNCand);
    
        PWM1T=patternNCandType(:)==(1:numAlph);
        PWMBgOut=sum(PWM1T);
        PWMBgOut=PWMBgOut./sum(PWMBgOut);
        
        PWMBg=PWMBg+alphab*(PWMBg0-PWMBgOut);
        PWMBg=PWMBg/sum(PWMBg);

    
        errorNormbg=norm(PWMBg0-PWMBgOut);
    
    
        errorNorm=errorNormbg;
    
        errorNormVec(iter, :)=errorNormbg;

        patternNCand=patternNCand(randperm(size(patternNCand, 1)), :);

        iter=iter+1;
    end
else
    errorNormVec=zeros(iterMax, 2);

    
    PWMBg0=PWMBg;
    
    while(errorNorm>0.001 && iter<iterMax)
        vertxType=embedPWM(patternNCand, PWMBg, vertxType, []);
        vertxType(patternCand(:))=patternCandTypes(:);

%         vertxType=embedPWM(patternCand, PWMEp, vertxType, patternCandTypes);
    
        zeroType=(vertxType==0);
        numZeros=sum(zeroType);
    
        if numZeros>0
    
            vertxType(zeroType)=randi(numAlph, 1, numZeros);
        end
    
    
    
    
        patternNCandType=vertxType(patternNCand);
    
        PWM1T=patternNCandType(:)==(1:numAlph);
        PWMBgOut=sum(PWM1T);
        PWMBgOut=PWMBgOut./sum(PWMBgOut);
    

    
        PWMBg=PWMBg+alphab*(PWMBg0-PWMBgOut);
        PWMBg=PWMBg/sum(PWMBg);

    
        errorNormbg=norm(PWMBg0-PWMBgOut);
    
    
        errorNorm=errorNormbg;

    
        errorNormVec(iter, :)=[0, errorNormbg];
        patternNCand=patternNCand(randperm(size(patternNCand, 1)), :);

        iter=iter+1;
    end
end

vertxType=vertxType(end-numVertex+1:end);


patternCandOut=patternCand+node0-1;

if options.isPlot
    figure
    plot(errorNormVec)
end
