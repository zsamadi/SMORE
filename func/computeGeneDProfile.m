function [pValMat, meDeltaMot]= computeGeneDProfile(geneEx, MotInTypes, pvalMin)


if nargin<3
    pvalMin=0.5;
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
    
    geneExDSi=geneExDS(:, iGS);

    probTwo=getLRCDF(geneExDSi,meDeltaMoti,nuMotifCell ,pvalMin);
    

    if abs(meDeltaMoti)>0 && probTwo>=0.5
        probTwoM=getLRCDF(geneExDSi,-meDeltaMoti,nuMotifCell ,pvalMin);
    else
        probTwoM=1;
    end






    pValMat(iGS)=min(probTwo, probTwoM);

end



