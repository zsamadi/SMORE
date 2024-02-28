function logPvalue=computePvalue(cntInfo, lenInfo, bernoulli)


if bernoulli<0

    nPosVec=repmat(lenInfo(1), length(cntInfo(:,1)), 1);

    nNegVec=repmat(lenInfo(2), length(cntInfo(:,1)), 1);


    logPvalue=getFisherLogPvalues(nPosVec, cntInfo(:,1), nNegVec, cntInfo(:,2));


elseif bernoulli>1

    meanLg=lenInfo(:, 1);
    stdLg=lenInfo(:, 2);
    stdLg(stdLg==0)=1;


    logPvalue= 1-logncdf(cntInfo(:, 1),meanLg,stdLg);
    logPvalue=max(logPvalue, eps);
    logPvalue=log(logPvalue);
    logPvalue=min(logPvalue, 0);

else

    cntP0F=cntInfo(:, 1)==0;

    pValAll=ones(size(cntInfo(:, 1)));

    cntInfoNZero=cntInfo(~cntP0F, :);

    pvalNZ = nbincdf(cntInfoNZero(:, 1)-1,cntInfoNZero(:, 2),1-bernoulli,'upper');

    pValAll(~cntP0F)=pvalNZ;


    logPvalue=log(pValAll);

    isInfinit=isinf(logPvalue);

    logPvalue(isInfinit)=getBinomLogPvalues(cntInfo(isInfinit,1), cntInfo(isInfinit,2), bernoulli);

end

