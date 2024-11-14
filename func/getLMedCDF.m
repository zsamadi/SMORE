function probi=getLMedCDF(geneExAlliSi, medMoti, nuMotifCell)

% probability of less than equal


numAllCell=size(geneExAlliSi, 1);

medianNum=floor(nuMotifCell/2);
[geneExAlliSiu, ~, ju]=unique(geneExAlliSi);

if length(geneExAlliSiu)==1
    probi=1;
elseif medMoti<geneExAlliSiu(1)
    probi=0;
else

    cnti=accumarray(ju, 1);

    % medMoti=medMot(iGS);

    % indexi=find(geneExAlliSiu-medMoti>=0, 1);
    % if ~any(geneExAlliSiu==medMoti)
    %     indexi=indexi-1;
    % end
    indexi=find(geneExAlliSiu-medMoti<=0,1, 'last');
        

    cntLeft=sum(cnti(1:indexi));

    cntRight=sum(cnti(indexi+1:end));

    
    imd=1;
    totalCases=nchoosekLog(numAllCell,nuMotifCell);
    probiv=zeros(nuMotifCell-cntLeft-cntRight+1, 1)-totalCases;

    if  cntLeft==0
        probi=0;
    elseif cntRight==0
        probi=1;
   elseif cntRight<=cntLeft
            probi=0.5;
    else
      probi=hygecdf(medianNum-1,numAllCell,cntLeft,nuMotifCell, 'upper');
      if probi==-1
         for medianNumi=max(nuMotifCell-cntLeft, 0):min(medianNum, cntRight)
                probiv(imd)=nchoosekLog(cntRight,medianNumi)+nchoosekLog(cntLeft,nuMotifCell-medianNumi);    
               imd=imd+1;
         end  
            probiv=probiv-totalCases;
            probi=sum(exp(probiv));
      end


        % special case where number of motif cells is even and the two middle
        % velues mean is smaller than motif median

        cntDownN=cntLeft-cnti(indexi);



        if (mod(nuMotifCell, 2)==0) && (cntDownN>(nuMotifCell-medianNum-2)) && (cntRight> (medianNum-1))

            if indexi>1
            if (geneExAlliSiu(indexi-1)+geneExAlliSiu(indexi))/2<medMoti && (geneExAlliSiu(indexi)+geneExAlliSiu(indexi+1))/2>medMoti
                probivSub=nchoosekLog(cntRight,medianNum)+nchoosekLog(cntDownN,nuMotifCell-medianNum-1)-totalCases;
                probi=probi-exp(probivSub);
            end
            end
        end


        probi=min(probi, 1);
    end
end