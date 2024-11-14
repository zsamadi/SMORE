function [pdfMDCell, edgeMDCell]= computeGeneDCDF(geneEx, MotInTypes)

% 
% maxGeneEx=max(geneEx);
% 
% geneEx=ceil(geneEx./max(geneEx)*100)/100*maxGeneEx;

numGenes=size(geneEx, 2);


pdfMDCell=cell(numGenes,1);
edgeMDCell=cell(numGenes,1);


geneExMed=median(geneEx);

geneExDS=geneEx-geneExMed;
geneExDS=sort(geneExDS);

nuMotifCell=sum(MotInTypes);

for iGS=1:size(geneExDS, 2)
    
    geneExDSi=geneExDS(:, iGS);
    geneExDSiU=unique(round(geneExDSi, 2));
    % if length(geneExDSiU)>100
    % 
    %     geneExDSiU=geneExDSiU(1:floor(length(geneExDSiU)/100):end);
    % end

    cdfMD=ones(length(geneExDSiU), 1);
    iU=1;
    probiLeft=0;
    while(probiLeft<0.99)
        meDeltaMoti=geneExDSiU(iU);
        
        probiLeft=getLOMedCDF(geneExDSi, meDeltaMoti, nuMotifCell);
        cdfMD(iU)=probiLeft;
        iU=iU+1;
        

    end

    lastiU=min(length(cdfMD), iU);
    cdfMD=cdfMD(1:lastiU);

    pdft=cdfMD-[0;cdfMD(1:end-1)];

    pdft=max(pdft, 0);
    pdft=pdft/sum(pdft);



    pdfMDCell{iGS}=pdft;
    edgeMDCell{iGS}=geneExDSiU(1:lastiU);

end




