function unqMer=getPSNmers(seqData, cTypes,options)

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
if options.isEnrich
    posFlag=posPWMS>0;
else
    posFlag=posPWMS>=max(posPWMS);
end

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

        nodesW=pSSites(iSdF, :);
        nodes=nodesW.';
        nodes=nodes(:);

        % if length(nodes)>5e6
        %     warning('\n input data is large:%d in getPSNmers \n',length(nodes));
        % end


        if ~isempty(nodes)

         isAllFixedType=options.cPTypesInit(nodes)==(options.fixedTypes(:)).';

        if all(isAllFixedType(:))
            iSit=false(length(nodes), 1);
            iSit(1)=true;
            % [~, iNDU]=unique(nodes{ii});
            % unqMersFlags(iNDU)=true;
        else
            iSit=getZNICIM(nodes,  ones(wMax, 1));
        end            
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





