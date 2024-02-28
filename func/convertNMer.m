function  [XNmers, merEyson]=convertNMer(seq,seqLens,cTypes, options)

W=options.W;
allMers=options.allMers;
if options.isNotNMer
    numShifts=seqLens-W+1;
    seqLensC=cumsum(seqLens);
    numShiftsRN=repelem([0;seqLensC(1:end-1)], numShifts);


    totShifts=sum(numShifts);
    numShiftsRi=(1:totShifts).';

    seqLensRi=[1;1+numShiftsRi(cumsum(numShifts(1:end-1)))];

    seqLensRi=repelem(seqLensRi, numShifts);

    seqLensRi=numShiftsRi-seqLensRi;


    merEyson=repelem(seqLensRi.', W,1);
    merEyson=merEyson(:);

    tmp2=repelem((1:W), totShifts, 1);
    tmp2=tmp2.';
    tmp2=tmp2(:);

    merEyson=(tmp2+merEyson);
    numShiftsRN=repelem(numShiftsRN, W, 1);
    merEyson=merEyson+numShiftsRN;

    XNmer=seq(merEyson);

    XNmer=(reshape(XNmer, W, totShifts)).';

    if allMers
        XNmerErased=false(totShifts, 1);
    else
        XNmerErased=any(XNmer==0, 2);
    end
    XNmer=XNmer(~XNmerErased, :);
    XNmerMaxWidth=repelem(seqLens, numShifts)-seqLensRi;
    XNmerMaxWidth=XNmerMaxWidth(~XNmerErased);
    XNmers=cTypes(XNmer);

    if options.isOKSingleNC
        XNmerSZFlag=(sum(XNmers==0, 2)==1);
        if any(XNmerSZFlag)
            XNmers(XNmerSZFlag, :)=options.cTypesInit(XNmer(XNmerSZFlag, :));
        end
    end

    if options.isOLess
        [XNmers, ~]=sort(XNmers, 2);
    end
    XNmers=[XNmers, XNmer, XNmerMaxWidth];

else
    XNmers=cTypes(seq(:, 1:end-1));

    if options.isOKSingleNC
        XNmerSZFlag=(sum(XNmers==0, 2)==1);
        if any(XNmerSZFlag)
            XNmers(XNmerSZFlag, :)=options.cTypesInit(seq(XNmerSZFlag, 1:end-1));
        end
    end


    if options.isOLess
        [XNmers, ~]=sort(XNmers, 2);
    end

    XNmers=[XNmers, seq];

    merEyson=(1:size(XNmers, 1)*W).';

end



