function pvalue=computeLogPvalue(cntInfo, lenInfo, testMethod)

% computes significance by using either fisher exact test or binomial test
% input arguments:
%   - cntInfo    : vector of posCount and negCount
%   - lenInfo    : vector of Np and NC, posLavg, negLavg average length of sequences, and W 
%   - testMethod : 'fisher' or 'binomial'
%  Output arguments:  
%   - pvalue     : computed pvalue

posCount=cntInfo(:,1);
negCount=cntInfo(:,2);

nPos=lenInfo(1);
nNeg=lenInfo(2);

if strcmpi(testMethod, 'fisher')
    pvalue=hygecdf(posCount-1,nNeg+nPos,nNeg,posCount+negCount,'upper'); 

    pvalueC=max(pvalue, realmin("double"));

    pvalue=log(pvalueC);

    nPosVec=repmat(nPos, length(posCount), 1);

    nNegVec=repmat(nNeg, length(posCount), 1);

    logPvalueMex=getLogFETPvalues(1, nPosVec, posCount, nNegVec, negCount);

    

    check=1;



%         logpi=zeros(nNeg-posCounti,1);

%         for i2=0:nNeg-posCounti-1
% 
%             logpi(i2+1,1)=logComb(nNeg, posCounti+i2)+logComb(nPos, negCounti-i2);
%         end
% 
%         logpimin=min(logpi);
%         logpiSum=logpimin+log(sum(exp(logpi-logpimin)));
%         logpv(i1)=logpiSum-logComb(nPos+nNeg, posCounti+negCounti);


    
else
    posLavg=lenInfo(3);
    negLavg=lenInfo(4);
    W=lenInfo(5);
    Pb=nPos*(posLavg-W+1)/(nPos*(posLavg-W+1)+nNeg*(negLavg-W+1));
    pvalue=1- binocdf(posCount-1,posCount+negCount,Pb);  
    pvalue=max(pvalue, realmin("double"));
end
