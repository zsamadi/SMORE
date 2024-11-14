function probi=getMedCDF(geneExAlliSi, medMoti, nuMotifCell)

% probability of greater than equal

numAllCell=size(geneExAlliSi, 1);

medianNum=floor(nuMotifCell/2);
[geneExAlliSiu, ~, ju]=unique(geneExAlliSi);

if length(geneExAlliSiu)==1
    probi=1;
elseif medMoti>geneExAlliSiu(end)
    probi=0;
else

    cnti=accumarray(ju, 1);

    % medMoti=medMot(iGS);

    indexi=find(geneExAlliSiu-medMoti>=0, 1);

    cntLeft=sum(cnti(1:indexi-1));

    cntRight=sum(cnti(indexi:end));

    
    imd=1;


    if cntRight==0
        probi=0;
    elseif cntLeft==0
        probi=1;
    elseif cntRight>=cntLeft
        probi=0.5;
    else
        probi=hygecdf(medianNum-1,numAllCell,cntRight,nuMotifCell, 'upper');

        if probi==0
            check=1;
        end


        if probi<-1
            totalCases=nchoosekLog(numAllCell,nuMotifCell);
            probiv=zeros(medianNum+1, 1)-totalCases;    
           for medianNumi=max(nuMotifCell-cntRight, 0):min(medianNum, cntLeft)
                probiv(imd)=nchoosekLog(cntLeft,medianNumi)+nchoosekLog(cntRight,nuMotifCell-medianNumi);
                imd=imd+1;
           end  
            probiv=probiv-totalCases;
            probit=sum(exp(probiv));
            if abs(probi-probit)/probi>0.01
                check=1;
            end


        end

        % special case where number of motif cells is even and the two middle
        % velues mean is smaller than motif median

        cntDownN=cntRight-cnti(indexi);



        if mod(nuMotifCell, 2)==0 && cntDownN>nuMotifCell-medianNum-2 && cntLeft> medianNum-1

            if (geneExAlliSiu(indexi-1)+geneExAlliSiu(indexi))/2<medMoti && (geneExAlliSiu(indexi)+geneExAlliSiu(indexi+1))/2>medMoti
                totalCases=nchoosekLog(numAllCell,nuMotifCell);
                probivSub=nchoosekLog(cntLeft,medianNum)+nchoosekLog(cntDownN,nuMotifCell-medianNum-1)-totalCases;
                probi=probi-exp(probivSub);
            end
        end


        probi=min(probi, 1);
    end
end