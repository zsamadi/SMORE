function pLRCDF=getLRCDF(geneExDSi,meDeltaMoti,nuMotifCell ,pvalMin)  

    probiRightPositive=getMedCDF(geneExDSi, meDeltaMoti, nuMotifCell);  

        if probiRightPositive<pvalMin
            probiLeftNegative=getLMedCDF(geneExDSi, -meDeltaMoti, nuMotifCell);
        else
            probiLeftNegative=0.5;
        end

        pLRCDF=probiRightPositive+probiLeftNegative;

		if pLRCDF<exp(-700) %  to prevent some rare special cases
		   pLRCDF=0.5;
		end
    
