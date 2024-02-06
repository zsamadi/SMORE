function logPvalue=computePvalue(cntInfo, lenInfo, bernoulli)

% computes significance by using either fisher exact test or binomial test
% input arguments:
%   - cntInfo    : vector of posCount and negCount
%   - lenInfo    : vector of Np and NC, posLavg, negLavg average length of sequences, and W 
%   - testMethod : 'fisher' or 'binomial'
%  Output arguments:  
%   - pvalue     : computed pvalue

% [~, idx]=sort(cntInfo(:,1), 'descend');
% 
% cntInfo=cntInfo(idx, :);

% Used in evaluateInitialSeeds, nestedSeedEnrichment, scoreModelPssm



if any(any(isnan(lenInfo)))
    check=1;
end

if bernoulli<0

    nPosVec=repmat(lenInfo(1), length(cntInfo(:,1)), 1);

    nNegVec=repmat(lenInfo(2), length(cntInfo(:,1)), 1);


    logPvalue=getFisherLogPvalues(nPosVec, cntInfo(:,1), nNegVec, cntInfo(:,2));


elseif bernoulli>1

        meanLg=lenInfo(:, 1);
        stdLg=lenInfo(:, 2);
        stdLg(stdLg==0)=1;

        if any(abs(imag(stdLg))>0)
            check=1;
        end
        % stdLg(meanLg==0)=1;


        logPvalue= 1-logncdf(cntInfo(:, 1),meanLg,stdLg);
        logPvalue=max(logPvalue, eps);
        logPvalue=log(logPvalue);
        logPvalue=min(logPvalue, 0);

else
    % cntInfo=round(cntInfo);

    cntP0F=cntInfo(:, 1)==0;

    pValAll=ones(size(cntInfo(:, 1)));

    cntInfoNZero=cntInfo(~cntP0F, :);
    % cntInfoNZero=max(cntInfoNZero, 1);


    pvalNZ = nbincdf(cntInfoNZero(:, 1)-1,cntInfoNZero(:, 2),1-bernoulli,'upper');

    pValAll(~cntP0F)=pvalNZ;


    logPvalue=log(pValAll);



    % logPvalue=max(logPvalue, -750);

    isInfinit=isinf(logPvalue);

    logPvalue(isInfinit)=getBinomLogPvalues(cntInfo(isInfinit,1), cntInfo(isInfinit,2), bernoulli);

    % 

    
%     posCount=cntInfo(:,1);
%     negCount=cntInfo(:,2);
%     nPos=lenInfo(1);
%     nNeg=lenInfo(2);
%     posLavg=lenInfo(3);
%     negLavg=lenInfo(4);
%     W=lenInfo(5);
%     Pb=nPos*(posLavg-W+1)/(nPos*(posLavg-W+1)+nNeg*(negLavg-W+1));
%     logPvalue=1- binocdf(posCount-1,posCount+negCount,Pb);  
%     logPvalue=log(max(logPvalue, realmin("double")));
end
% pvalue=gather(pvalue);

