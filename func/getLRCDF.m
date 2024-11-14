function pLRCDF=getLRCDF(geneExDSi,meDeltaMoti,nuMotifCell ,pvalMin)  

    probiRightPositive=getMedCDF(geneExDSi, meDeltaMoti, nuMotifCell);  

        if probiRightPositive<pvalMin
            probiLeftNegative=getLMedCDF(geneExDSi, -meDeltaMoti, nuMotifCell);
        else
            probiLeftNegative=0.5;
        end

        pLRCDF=probiRightPositive+probiLeftNegative;



    pLRCDF=max(pLRCDF, exp(-700));
    