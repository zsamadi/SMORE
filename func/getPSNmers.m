function unqMer=getPSNmers(seqData, cTypes,options)


% Counts all WMers in the input primary and control data
% input arguments:
%   - seqData:
%       - seqData.pSeq  : primary sequence
%       - seqData.nSeq  : control sequence
%       - seqData.pHSeq : primary holdoput sequence
%       - seqData.nHSeq : control holdoput sequence
%       - seqData.PWM0  : backgroung model
%       - seqData.lens  : data lengths [size(posSeq90, 1),size(negSeq90, 1),size(pHoldSeq, 1),size(nHoldSeq, 1)];
%   - options
%       - options.NEVAL     : Number of motifs for initial evaluation
%       - options.NREF      : Number of motifs for refinement
%       - options.prior     : Dirichlet prior for initial evaluation
%       - options.REFPrior  : Dirichlet prior for refinement step
%  Output arguments:
%   - infoCell
%     - uniMers     : all seeds of length W, counts of Wmer seeds in primary and control data
%     - primCntsF   : counts of Wmer seeds in primary data
%     - ctrlCntsF   : counts of Wmer seeds in control data
%   - wmerZoops     : structure of cells containing zoops WMers of primary and control
%       - wmerZoops.pos=posWmerZoopsCell;
%       - wmerZoops.neg=negWmerZoopsCell;
% wMin=options.wMin;
wMax=options.wMax;
numCells=options.numCells;


cnvOptions.W=wMax;
cnvOptions.allMers=false;
cnvOptions.isOLess=options.isOLess;
cnvOptions.isOKSingleNC=options.isOKSingleNC;


cTypesEx=[cTypes(:); numCells+1;0];



if options.isPos


    pWeights=seqData.posWeight;
    cnvOptions.cTypesInit=[options.cPTypesInit; numCells+1;0];






else


    pWeights=seqData.negWeight;
    cnvOptions.cTypesInit=[options.cNTypesInit; numCells+1;0];

end
% Sorts input sequence matrix in NMers and Zoops NMers



cnvOptions.isNotNMer=false;
pSSites=convertNMer(seqData.pSeq,[],cTypesEx, cnvOptions);

sites=[(1:size(pSSites, 1)/2).';(1:size(pSSites, 1)/2).'];
sites=sites(:);

[~, iStW]=sortrows([pSSites(:, 1:wMax),-pWeights,sites, -pSSites(:, end)]);


pSSites=pSSites(iStW, :);

pWeights=pWeights(iStW);



PWMS=options.PWMS;

[pSSitesU, ~, jPSU]=unique(pSSites(:, 1:wMax), 'rows');


posPWMS=scoreWords(pSSitesU  ,PWMS, options);

posPWMS=posPWMS(jPSU);

% posPWMSv=scoreWords(pSSites(:, 1:wMax)  ,PWMS(:, end:-1:1), options);
% posPWMS=max(posPWMSt, posPWMSv);

posFlag=posPWMS>0;

if any(posFlag)

    posPWMS=posPWMS(posFlag);

    poSeeds=pSSites(posFlag, 1:wMax);

    [poSeedsU, iU, jU]=unique(poSeeds, 'rows');
    posPWMS=posPWMS(iU);

    pSSites=pSSites(posFlag, wMax+1:2*wMax);

    pSWeights=pWeights(posFlag);

    pSSiteCell=cell(size(poSeedsU, 1), 1);
    pSWeightCell=cell(size(poSeedsU, 1), 1);

    for iSd=1:size(poSeedsU, 1)
        iSdF=(jU==iSd);
        weighsAll=pSWeights(iSdF, :);
        % [weighsAll, ist]=sort(weighsAll, 'descend');

        nodesW=pSSites(iSdF, :);
        % nodesW=nodesW(ist, :);
        nodes=nodesW.';
        nodes=nodes(:);

        if length(nodes)>1e6
            warning('come on lad, input data is too big:%d in getPSNmers \n',length(nodes));
        end


        if ~isempty(nodes)
            iSit=getZNICIM(nodes,  ones(wMax, 1));
            iSit=sum(reshape(iSit, wMax, []));
            iSit=iSit.';
            iSit=iSit>0;
    
            pSSiteCell{iSd}=nodesW(iSit, :);
    
            pSWeightCell{iSd}=weighsAll(iSit, :);
        else
            pSSiteCell{iSd}=[];

            pSWeightCell{iSd}=[];
        end

    end

    [posPWMS, ist]=sort(posPWMS, 'descend');


    unqMer.posPWMS=posPWMS;

    unqMer.seeds=poSeedsU(ist, :);

    unqMer.siteCell=pSSiteCell(ist);

    unqMer.weightCell=pSWeightCell(ist);
else
    unqMer.posPWMS=[];

    unqMer.seeds=[];

    unqMer.siteCell={[]};

    unqMer.weightCell={[]};
end





