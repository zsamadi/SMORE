function  [XNmers, merEyson]=convertNMer(seq,seqLens,cTypes, options)

% Sorts input sequence matrix in NMers and Zoops NMers
% Input arguments:
%   - seq   : input matrix of sequences
%   - W     : length of mers
%   - zpsOut: true or false, if true, outputs zoops nmenrs as well
% Output arguments:
%   - XZoops: Zoops WMers, 0 if zpsOut is false
%   - XNmers: output WMers
%used in countSeed and scoreSeq and seqFilterNew


% numPSeqs=options.numPSeqs;
W=options.W;
% rvPath=options.rvp;
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


    % numShiftsRN=numShiftsRN;
    numShiftsRN=repelem(numShiftsRN, W, 1);
    merEyson=merEyson+numShiftsRN;
    % XNmer=zeros(W, totShifts);
    % for ii=1:W
    %     XNmer(ii, :)=seq(merEyson(ii:W:end));
    % end
    % XNmer=XNmer.';

    XNmer=seq(merEyson);

    XNmer=(reshape(XNmer, W, totShifts)).';

    if allMers
        XNmerErased=false(totShifts, 1);
    else
        XNmerErased=any(XNmer==0, 2);
    end
    XNmer=XNmer(~XNmerErased, :);

    % if (rvPath)
    % %     XNmerSitesP=repelem((1:numPSeqs).', numShifts(1:numPSeqs), 1);
    % %     XNmerSitesP=XNmerSitesP(:);
    %     if options.pnSeq
    % %         numNSeqs=N/2-numPSeqs;
    %
    % %         XNmerSitesN=repelem((numPSeqs+1:N/2).', numShifts(2*numPSeqs+1:2*numPSeqs+numNSeqs), 1);
    % %         XNmerSitesN=XNmerSitesN(:);
    % %         XNmerSites=[XNmerSitesP;XNmerSitesP;XNmerSitesN;XNmerSitesN];
    %
    %     else
    % %         XNmerSites=[XNmerSitesP;XNmerSitesP];
    %     end
    %
    %
    % else
    %
    % %     XNmerSites=repelem((1:N).', numShifts, 1);
    %
    % %     XNmerSites=XNmerSites(:);
    % end



    XNmerMaxWidth=repelem(seqLens, numShifts)-seqLensRi;
    XNmerMaxWidth=XNmerMaxWidth(~XNmerErased);



    % XNmerSites=XNmerSites(~XNmerErased);


    % XNmerV=XNmer(:);

    % XNmers=XNmerV;


    % XNmers(XNmerV>0)=cTypes(XNmerV(XNmerV>0));
    XNmers=cTypes(XNmer);

    if options.isOKSingleNC
        XNmerSZFlag=(sum(XNmers==0, 2)==1);
        if any(XNmerSZFlag)
            XNmers(XNmerSZFlag, :)=options.cTypesInit(XNmer(XNmerSZFlag, :));
        end
    end



    % XNmers=reshape(XNmers, size(XNmer));

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



    % XNmers=reshape(XNmers, size(XNmer));

    if options.isOLess
        [XNmers, ~]=sort(XNmers, 2);
    end



    XNmers=[XNmers, seq];

    merEyson=(1:size(XNmers, 1)*W).';




end



