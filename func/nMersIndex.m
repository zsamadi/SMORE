function merEyson=nMersIndex(seqLens, W)
numShifts=seqLens-W+1;

seqLensC=cumsum(seqLens);
numShiftsRN=repelem([0;seqLensC(1:end-1)], numShifts,1);


totShifts=sum(numShifts);
numShiftsRi=(1:totShifts).';

numShiftsRit=[1;1+numShiftsRi(cumsum(numShifts(1:end-1)))];

numShiftsR=repelem(numShiftsRit, numShifts,1);
seqLensRi=numShiftsRi-numShiftsR;

tmp=repelem(seqLensRi.', W,1);
tmp=tmp(:);

tmp2=repelem((1:W), totShifts, 1);
tmp2=tmp2.';
tmp2=tmp2(:);

merEyso=(tmp2+tmp);


sitest=numShiftsRN;
sitesE=repelem(sitest, W, 1);
merEyson=merEyso+sitesE;