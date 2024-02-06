function [unqMers, nMersMaxOut, unqMersFlagsOut]=countWSeeds(nMers, W,nMersMax,unqMersFlagsIn, options)
numCells=options.numCells;
wMax=options.wMax;
numNodes=options.numNodes;
pnWeights=options.pnWeights;

% tmp=nMers;


% sites=(1:size(nMers, 1));


if options.rvp
    if options.isPN
        numSeqScale=(wMax-options.MinSeedWidth+1)*options.numSeqs;
        sitesP=[(1:numSeqScale).';(1:numSeqScale).'];

        sitesNL=size(nMers, 1)-size(sitesP, 1);


        sites=[sitesP;(1:sitesNL/2).';(1:sitesNL/2).'];
    else
        sites=[(1:size(nMers, 1)/2).';(1:size(nMers, 1)/2).'];
    end        
end
sites=sites(:);

% nMersW=[nMers(:, 1:W), nMers(:, end-wMax:end-wMax+W-1), nMers(:, end)];
% 
% nMersW(:, end)=-nMersW(:, end);
% [~, iStW]=sortrows([nMers(:, 1:W), nMers(:, end-wMax:end-wMax+W-1), -nMers(:, end)]);
% [~, iStW]=sortrows([nMers(:, 1:W),pnWeights,sites, -nMers(:, end)]);
[~, iStW]=sortrows([nMers(:, 1:W),-pnWeights,sites, -nMers(:, end)]);

nMers=nMers(iStW, :);

pnWeights=pnWeights(iStW);

isPlot71Contacts=false;

if isPlot71Contacts
    aat=nMers(:, 1:2);
    aatf=any(aat==71, 2);
    aan=nMers(:, 3:4);
    
    pnWeights71=pnWeights(aatf);
    aat71=aat(aatf, :);
    
    aan71=aan(aatf, :);
    
    aat71=sort(aat71, 2);
    
    [aat71U, iU, jU]=unique(aat71, 'rows');
    
    aaw71=accumarray(jU, pnWeights71);
    
    figure
    plot(aaw71/sum(aaw71)*100)
    hold on
end




nMersAndSitesSort=[nMers(:, 1:W), nMers(:, end-wMax:end-wMax+W-1)];

% [nMersAndSitesSort, iZSort]=sortrows(nMersAndSitesSort);

% nMers=nMers(iZSort, :);

%%% C IMPLEMENTATION

nodes=nMersAndSitesSort(:, W+1:end);
nodes=nodes.';
nodes=nodes(:);
unqEMers=nMersAndSitesSort(:, 1:W);
nodesFlag=sum(unqEMers.*(numCells+2).^(W-1:-1:0), 2);
if ~isempty(unqMersFlagsIn)
    [nodesFlagU, iNodesU]=unique(nodesFlag);
else
    [~, iNodesU]=unique(nodesFlag);
end

% pWeights=pnWeights;

% pWeights(nMersAndSitesSort(:, W+1)>numNodes)=0;
% pCWeights=accumarray(jNodesU, pWeights>0);
% pWeights=accumarray(jNodesU, pWeights);
% pCWeights(pCWeights==0)=1;
% pWeights=[pWeights./pCWeights,pCWeights];


% nWeights=pnWeights;

% nWeights(nMersAndSitesSort(:, W+1)<=numNodes)=0;
% nCWeights=accumarray(jNodesU, nWeights>0);
% nWeights=accumarray(jNodesU, nWeights);
% nCWeights(nCWeights==0)=1;
% nWeights=[nWeights./nCWeights, nCWeights];





unqEMers=unqEMers(iNodesU,:);
nMersMaxOut=nMers(iNodesU, :);

iNodesU0=[iNodesU(2:end);length(nodesFlag)+1];
iNodesUD=iNodesU0-[1; iNodesU0(1:end-1)];

pnWeights=mat2cell(pnWeights, iNodesUD,1);




iNodesU=(iNodesU-1)*W+1;

iNodesU=[iNodesU(2:end);W*length(nodesFlag)+1];

iNodesUD=iNodesU-[1; iNodesU(1:end-1)];
nodes=mat2cell(nodes, iNodesUD,1);


uniqueMersV=(unqEMers(:, end)<numCells+1);
uniqueMersV0=~any(unqEMers==0, 2);
uniqueMersV=uniqueMersV&uniqueMersV0;
unqEMers=unqEMers(uniqueMersV, :);
% pWeights=pWeights(uniqueMersV, :);
% nWeights=nWeights(uniqueMersV, :);

pnWeights=pnWeights(uniqueMersV);


nodes=nodes(uniqueMersV);
nNodesU=sum(uniqueMersV);


pCounts=zeros(nNodesU, 1);
nCounts=zeros(nNodesU, 1);

pWeights=zeros(nNodesU, 1);
nWeights=zeros(nNodesU, 1);

iivd=zeros(numNodes,1);
iid=1;
iU=zeros(nNodesU,1);
sitesCell=cell(nNodesU,1);

weightsCell=cell(nNodesU,1);

unqMersFlagsC=cell(nNodesU, 1);

if ~isempty(unqMersFlagsIn)

    nodesFlagU=nodesFlagU(uniqueMersV);
end

% nMersW=nMersAndSitesSort(:, 1:W);

MersV=(nMersAndSitesSort(:, W)<numCells+1);
MersV0=~any(nMersAndSitesSort(:, 1:W)==0, 2);
MersV=MersV&MersV0;


nMersAndSitesSort=nMersAndSitesSort(MersV, :);

% indexCheck=find(all(sort(unqEMers, 2)==[44, 71], 2));


if ~isempty(unqMersFlagsIn)
    [~, locb]=ismember(nodesFlagU, unqMersFlagsIn.nodesFlagU);

    for ii=1:nNodesU



        cNodes=nodes{ii};

        locbii=locb(ii);


        if locb(ii)
            cNodesP=unqMersFlagsIn.nodesC{locbii};
        else
            cNodesP=[];
        end

        if length(cNodes)==length(cNodesP)
            if sum(abs(cNodes-cNodesP))==0
                unqMersFlagsC{ii}=unqMersFlagsIn.unqMersFlagsC{locbii};
                unqMersFlagsC{ii}=unqMersFlagsIn.unqMersFlagsC{locbii};
                cNodesW=reshape(cNodes, W, []);
                cNodesW=cNodesW.';

                cNodes0=cNodesW(:,1);
                cNodes0=cNodes0(unqMersFlagsC{ii}>0);
                pCounts(ii)=sum(cNodes0<=numNodes);
                nCounts(ii)=sum(cNodes0>numNodes);
                iU(ii)=sum(unqMersFlagsC{ii});
                sitesCell{ii}=cNodesW(unqMersFlagsC{ii}>0, :);
            else

                iivd(iid)=ii;
                iid=iid+1;

            end

        else
            iivd(iid)=ii;
            iid=iid+1;
        end


    end

    iivd=iivd(1:iid-1);


    for iidn=1:length(iivd)

        ii=iivd(iidn);

        cNodes=nodes{ii};
        


        unqMersFlags=getZNICIM(cNodes,  ones(W, 1));
        unqMersFlags=reshape(unqMersFlags, W, []);
        unqMersFlags=sum(unqMersFlags);
        unqMersFlags=unqMersFlags.';
        unqMersFlagsC{ii}=unqMersFlags;
        cNodesW=reshape(cNodes, W, []);
        cNodesW=cNodesW.';
        cNodes0=cNodesW(:, 1);
        cNodes0=cNodes0(unqMersFlags>0);
        pCounts(ii)=sum(cNodes0<=numNodes);
        nCounts(ii)=sum(cNodes0>numNodes);

        iU(ii)=sum(unqMersFlags);
        sitesCell{ii}=cNodesW(unqMersFlags>0, :);



    end

    jU=(1:length(iU)).';

    jU=repelem(jU, iU,1);

    iU=cumsum([1;iU(1:end-1)]);

    unqMersFlagsAll=vertcat(unqMersFlagsC{:});



    %     nMerZoops=nMersAndSitesSort(unqMersFlagsAll>0, :);

    % iZ=iZSort(unqMersFlagsAll>0);




else

    siteNum0=1;
    for ii=1:nNodesU
        % cNodes=nodes{ii};
        % 
        % if any(indexCheck==ii)
        %     check=1;
        % end

        unqMersFlags=getZNICIM(nodes{ii},  ones(W, 1));

        % just for testing to see what happens if we don't consider znics

        % unqMersFlags=unqMersFlags+1;
        unqMersFlags=reshape(unqMersFlags, W, []);
        unqMersFlags=sum(unqMersFlags);
        unqMersFlags=unqMersFlags.';
        unqMersFlags=unqMersFlags>0;


        unqMersFlagsC{ii}=unqMersFlags;
        cNodesW=reshape(nodes{ii}, W, []);
        cNodesW=cNodesW.';
        cNodesW=cNodesW(unqMersFlags, :);
        % cNodes0=cNodesW(:, 1);
        % cNodes0=cNodes0(unqMersFlagsC{ii}>0);
        pCounts(ii)=sum(cNodesW(:, 1)<=numNodes);
        nCounts(ii)=sum(cNodesW(:, 1)>numNodes);
        
        % sitesCell{ii}=cNodesW(unqMersFlags>0, :);

        pnWeightsi=pnWeights{ii};
        pnWeightsi=pnWeightsi(unqMersFlags);

        pWeights(ii)=sum(pnWeightsi(cNodesW(:, 1)<=numNodes));
        nWeights(ii)=sum(pnWeightsi(cNodesW(:, 1)>numNodes));



    end
    unqMersFlagsAll=vertcat(unqMersFlagsC{:});

    nMersAndSitesSort=nMersAndSitesSort(unqMersFlagsAll, :);

    nMersAndSitesSort=nMersAndSitesSort(:, 1:W);

    if ~isempty(nMersAndSitesSort)
        [~,iU, jU]=uniqueSorted(nMersAndSitesSort, numCells+2);
    else
        iU=[];
        jU=[];
    end





end


if isPlot71Contacts
    
    aaUM=unqEMers;
    aaUMf=any(aaUM==71, 2);
    
    pnUMWeights71=pWeights(aaUMf);
    aaUM71=aaUM(aaUMf, :);
    
    aaUM71=sort(aaUM71, 2);
    
    [aaUM71U, iUM, jUM]=unique(aaUM71, 'rows');
    
    aawUM71=accumarray(jUM, pnUMWeights71);
    
    
    
    plot(aawUM71/sum(aawUM71)*100)
    grid on
    xlabel('cell type')
    ylabel('connection percentage')
    legend('raw weights', 'unique weights')
end




% unqMersFlagsOut.unqMersFlagsC=unqMersFlagsC;
% unqMersFlagsOut.nodesFlagU=nodesFlagU;
% unqMersFlagsOut.nodesC=nodes;

unqMersFlagsOut=[];

XNmerzMaxWidth=nMers(:,end);
XNmerzMaxWidth=XNmerzMaxWidth(unqMersFlagsAll>0);





% jUP=jU(nMerZoops(:,W+1)<=numNodes);
% jUN=jU(nMerZoops(:,W+1)>numNodes);
% nuMers=size(unqEMers, 1);
% pCounts=zeros(nuMers, 1);
% nCounts=zeros(nuMers, 1);

% if ~isempty(jUP)
%     pCountst = accumarray(jUP,1);
%     pCounts(1:length(pCountst))=pCountst;
%
% end
% if ~isempty(jUN)
%
%     nCountst = accumarray(jUN,1);
%     nCounts(1:length(nCountst))=nCountst;
%
% end





%
% pCounts=pCounts(uniqueMersV);
% nCounts=nCounts(uniqueMersV);
%
% % BE CAREFULL
% nCounts=round(nCounts/sum(nCounts)*sum(pCounts));

unqMers.counts=[pCounts, nCounts];

unqMers.weights=[pWeights, nWeights];


mWDTSort=sortrows([jU, XNmerzMaxWidth]);
if ~isempty(mWDTSort)
    XNmerzWidth=mWDTSort(:,2);
    XNmerzMaxWidth=XNmerzWidth(iU);
else
    XNmerzMaxWidth=[];
end

% XNmerzMaxWidth=XNmerzMaxWidth(uniqueMersV);
% STREME goes to this depth
% unqMers.maxWidth=min(XNmerzMaxWidth, wMax+1);
unqMers.maxWidth=XNmerzMaxWidth;


% siteCell is later used in initial and nested seed enrichment, and modify counts
% sitesAll=nMerZoops(:, W+1:end);

% iUD=[iU(2:end);size(nMerZoops,1)+1]-iU;
% this will be used in case of modify count
% sitesCell=mat2cell([sitesAll, XNmerzWidth],iUD);

% sitesCell=mat2cell(sitesAll, iUD,W);
% for iiD=1:length(iUD)
%     sitesi=sitesAll(iiD0:iiD0+iUD(iiD)-1, :);
%
%
%     sitesCell{iiD}=sitesi;
%     iiD0=iiD0+iUD(iiD);
% end

% sitesCell=sitesCell(uniqueMersV);


% unqMers.siteCell=sitesCell;

%
% sitesCell=sitesCellAll(uniqueMersV);
% unqMers.siteCell=sitesCell;
unqMers.w=W;



if (W<wMax && (~isempty(nMersMax)))

    nMersMaxV=nMersMax(:, W)<(numCells+1);
    nMersMaxV0=~any(nMersMax(:, 1:W)==0, 2);
    nMersMaxV=nMersMaxV&nMersMaxV0;
    nMersMax=nMersMax(nMersMaxV, :);


    [nMersMaxT, iSt]=sortrows([nMersMax(:, 1:W),nMersMax(:, end-wMax:end-wMax+W-1), -nMersMax(:, end)]);

    [~,iT]=uniqueSorted(nMersMaxT(:, 1:W), numCells+1);
    unqMers.stExt=nMersMax(iSt, 1:end-1);

    unqMers.seeds=unqMers.stExt(:, 1:wMax);

    iTN=min(iT+1, [iT(2:end)-1;length(iSt)]);

    unqMers.seeds=unqMers.seeds(iTN, :);
    unqMers.stExt=unqMers.stExt(iTN,wMax+1:end);

    % if size(unqEMers, 1)~=size(unqMers.seeds, 1)
    %     check=1;
    % end

else

    unqMers.seeds=unqEMers;
    unqMers.seedsIndex=sum(unqEMers.*(numCells.^(W-1:-1:0)), 2);


end


% check=1;

