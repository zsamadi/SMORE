function [pValMat, meDeltaMot]= computeGeneDProfile(geneEx, MotInTypes, isTwoSided)


if nargin<3
    isTwoSided=true;
end



numGenes=size(geneEx, 2);


pValMat=zeros(numGenes,1);



geneExMot=geneEx(MotInTypes, :);
geneExMed=median(geneEx);

geneExDS=geneEx-geneExMed;
geneExDS=sort(geneExDS);
meDeltaMot=median(geneExMot)-geneExMed;

nuMotifCell=sum(MotInTypes);

for iGS=1:size(geneExDS, 2)


    meDeltaMoti=meDeltaMot(iGS);
    if meDeltaMoti~=0 && mod(nuMotifCell, 2)==0
        check=1;
    end

    if iGS==155
        check=1;
    end
    
    geneExDSi=geneExDS(:, iGS);
    probiRightPositive=getMedCDF(geneExDSi, meDeltaMoti, nuMotifCell);


    probiLeftNegative=getLMedCDF(geneExDSi, -meDeltaMoti, nuMotifCell);

    if isTwoSided
        probTwo=probiRightPositive+probiLeftNegative;
    else
        probTwo=probiRightPositive;
    end


    probiLeftPositive=getLMedCDF(geneExDSi, meDeltaMoti, nuMotifCell);


    probiRightNegative=getMedCDF(geneExDSi, -meDeltaMoti, nuMotifCell);
    % probiLeftPrime=getMedCDF(geneExDSi, -meDeltaMoti, nuMotifCell);

    if isTwoSided
        probTwoM=probiRightNegative+probiLeftPositive;
    else
        probTwoM=probiLeftPositive;
    end

    if abs(log(min(probTwo, probTwoM)))>20
        check=1;
    end



    pValMat(iGS)=min(probTwo, probTwoM);

     if abs(log(pValMat(iGS)))>20
         check=1;
     end

end

check=1;



