function pValMat=createGeneEProfile(geneEx, MotInTypes)

isReplace=true;
geneExMot=geneEx(MotInTypes, :);
geneExDMot=geneEx;


iterMax=1000;
numBins=100;
medType=median(geneExDMot);

meDeltaMot=median(geneExMot)-medType;



numCells=size(geneExMot, 1);

numGenes=size(geneEx, 2);
medianMat=zeros(iterMax, numGenes);



for iterProf=1:iterMax

    randMat=rand(size(geneEx));

    [~, randInd]=sort(randMat);

    randInd=randInd(1:numCells, :);



    
    allIDX=(1:numel(geneEx));
    % randInd=randi(size(geneEx, 1), numCells, size(geneEx, 2));
    randInd=randInd+(0:size(geneEx, 1):numel(geneEx)-1);

    randInd=randInd(:);


    geneExi=geneEx(:);
    geneExi=geneExi(randInd);
    geneExi=reshape(geneExi, [], numGenes);

    if ~isReplace
        allIDX(randInd)=[];
        geneExDi=geneEx(:);
        geneExDi=geneExDi(allIDX);
        geneExDi=reshape(geneExDi, [], numGenes);
        dmediani=median(geneExi)-median(geneExDi);
    else
        dmediani=median(geneExi);

    end

    medianMat(iterProf, :)=dmediani;

end

if isReplace
    medianMat=medianMat-medType;
end



pValMat=zeros(numGenes,1);



for iGene=1:numGenes

    
  
    [NMed,edgeMed] = histcounts(medianMat(:, iGene), numBins);



    pvalPositive=getDeltaPvalue(NMed,edgeMed, meDeltaMot(iGene));



    pvalNegetive=2-pvalPositive;



    pValMat(iGene)=min(pvalPositive, pvalNegetive);





end





