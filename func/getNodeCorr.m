function [numMNodes,numPMNodes]=getNodeCorr(ppHMNodesCell,extMotif,cPTypes,nodeSID, fixedTypes,W,nMCorrv)
nodeSIDT=unique(nodeSID);
numSIDT=length(nodeSIDT);
nMCorr=length(nMCorrv);
numMNodes=zeros(nMCorr, numSIDT);
numPMNodes=zeros(nMCorr, numSIDT);
cPTypes=cPTypes(:);
motifPvalV=zeros(nMCorr, 1);

for imti=1:nMCorr
    imt=nMCorrv(imti);
    mNodesi=ppHMNodesCell{imt};
   %  cNodes=mNodesi.';
   %  cNodes=cNodes(:);
   % 
   %  fixedNodesFlag=cPTypes(cNodes)==(fixedTypes(:)).';
   %  fixedNodesFlag=fixedNodesFlag(:);
   % 
   % 
   %  if sum(fixedNodesFlag)/length(fixedNodesFlag)<0.25
   % 
   % 
   % unqMersFlags=getZNICIM(cNodes,  ones(W, 1));
   %      unqMersFlags=reshape(unqMersFlags, W, []);
   %      unqMersFlags=sum(unqMersFlags);
   %      unqMersFlags=unqMersFlags.';
   %      cNodesW=reshape(cNodes, W, []);
   %      cNodesW=cNodesW.';
   %      mNodesi=cNodesW(:, 1);
   %      mNodesi=mNodesi(unqMersFlags>0);

    mNodesSIDi=nodeSID(mNodesi(:, 1));
    [mNodesSIDiU, ~, jSIDU]=unique(mNodesSIDi);
    cSIDi=accumarray(jSIDU, 1);
    % cSIDi(cSIDi<5)=0;
    exMotifi=extMotif{imt};
    % 
    % cntiT(mNodesSIDiU)=cSIDi*abs(exMotifi.testPvalue);
    UId=ismember(nodeSIDT, mNodesSIDiU);

    numMNodes(imti,UId) =cSIDi;
    % numPMNodes(imt, UId)=cSIDi*abs(exMotifi.testPvalue);
    numPMNodes(imti, UId)=cSIDi;

    motifPvalV(imti)=abs(exMotifi.testPvalue);

    % end
end



% numMNodes=numMNodes./sum(numMNodes);
% numPMNodes=numPMNodes./sum(numPMNodes);



% numMNodes(:, 1:end-1)=numMNodes(:, 1:end-1)./sum(numMNodes(:, 1:end-1), 2);
% numPMNodes(:, 1:end-1)=numPMNodes(:, 1:end-1)./sum(numPMNodes(:, 1:end-1), 2);

numMNodes=numMNodes./sum(numMNodes, 2);
% numPMNodes=numPMNodes./sum(numPMNodes, 2)./nMCorrv(:);
numPMNodes=numPMNodes./sum(numPMNodes, 2).*motifPvalV;

% numMNodes=numMNodes./sum(numMNodes);
% numPMNodes=numPMNodes./sum(numPMNodes);


numMNodes=numMNodes(~any(isnan(numMNodes), 2), :);
numPMNodes=numPMNodes(~any(isnan(numPMNodes), 2), :);

% 
% numMNodes0L=numMNodes;
% 
% numMNodes0L(numMNodes0L==0)=1;
% 
% etr0py=-sum(numMNodes(:, 1:end-1).*log2(numMNodes0L(:, 1:end-1)), 2);
% 
% numPMNodes0L=numPMNodes;
% 
% numPMNodes0L(numPMNodes0L==0)=1;
% 
% etr1py=-sum(numPMNodes(:, 1:end-1).*log2(numPMNodes0L(:, 1:end-1)), 2);
% 
% 
% [~, numMNodeSi]=sortrows([numMNodes(:, end), etr0py], 'descend');
% 
% 
% numMNodes=numMNodes(numMNodeSi, :);
% etr0py=etr0py(numMNodeSi);
% 
% outSize=max(find(etr0py==0, 1, 'first')-1, 2*find(numMNodes(:, end)>0, 1, 'last'));
% % outSize=find(numMNodes(:, end)>0, 1, 'last');
% if ~isempty(outSize)
%     outSize=min(outSize, size(numMNodes, 1));
% else
%     outSize=size(numMNodes, 1);
% end
% 
% numMNodes=numMNodes(1:outSize, :);
% 
% 
% [~, numPMNodeSi]=sortrows([numPMNodes(:, end), etr1py], 'descend');
% 
% numPMNodes=numPMNodes(numPMNodeSi, :);
% 
% numPMNodes=numPMNodes(1:outSize, :);
