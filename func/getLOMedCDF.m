function probi=getLOMedCDF(geneExAlliSi, medMoti, nuMotifCell)

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

    indexi=find(geneExAlliSiu-medMoti<=0,1, 'last');
        

    cntLeft=sum(cnti(1:indexi));

    cntRight=sum(cnti(indexi+1:end));

    probiv=zeros(medianNum+1, 1);
    imd=1;
    totalCases=nchoosekLog(numAllCell,nuMotifCell);

    if  cntLeft==0
        probi=0;
    elseif cntRight==0
        probi=1;
    else

      probi=hygecdf(medianNum,numAllCell,cntLeft,nuMotifCell, "upper");
      % if probi==0
      % 
      %    for medianNumi=max(nuMotifCell-cntLeft, 0):min(medianNum, cntRight)
      %           probiv(imd)=nchoosekLog(cntRight,medianNumi)+nchoosekLog(cntLeft,nuMotifCell-medianNumi);    
      %          imd=imd+1;
      %    end  
      %       probiv=probiv-totalCases;
      %       probi=sum(exp(probiv));
      % end


        % special case where number of motif cells is even and the two middle
        % velues mean is smaller than motif median

        cntDownN=cntLeft-cnti(indexi);



        if (mod(nuMotifCell, 2)==0) && (cntDownN>(nuMotifCell-medianNum-2)) && (cntRight> (medianNum-1))

            if indexi>1
            if (geneExAlliSiu(indexi-1)+geneExAlliSiu(indexi))/2<medMoti
                probivSub=nchoosekLog(cntRight,medianNum)+nchoosekLog(cntDownN,nuMotifCell-medianNum-1)-totalCases;
                probi=probi-exp(probivSub);
            end
            end
        end


        probi=min(probi, 1);
    end
end