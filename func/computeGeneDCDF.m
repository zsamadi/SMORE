function [pdfMDCell, edgeMDCell]= computeGeneDCDF(geneEx, MotInTypes)


numGenes=size(geneEx, 2);


pdfMDCell=cell(numGenes,1);
edgeMDCell=cell(numGenes,1);


geneExMed=median(geneEx);

geneExDS=geneEx-geneExMed;
geneExDS=sort(geneExDS);

nuMotifCell=sum(MotInTypes);

for iGS=1:size(geneExDS, 2)
    
    geneExDSi=geneExDS(:, iGS);
    geneExDSiU=unique(geneExDSi);
    cdfMD=zeros(length(geneExDSiU), 1);

    for iU=1:length(geneExDSiU)
        meDeltaMoti=geneExDSiU(iU);
        probiLeft=getLMedCDF(geneExDSi, meDeltaMoti, nuMotifCell);
        cdfMD(iU)=probiLeft;

    end

    pdfMDCell{iGS}=cdfMD-[0;cdfMD(1:end-1)];
    edgeMDCell{iGS}=geneExDSiU;

end




