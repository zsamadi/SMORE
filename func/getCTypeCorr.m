
function numCTNodes=getCTypeCorr(cPTypes,nodeSID, fixedTypes)

nodeSIDT=unique(nodeSID);
numSIDT=length(nodeSIDT);

cPTypesU=unique(cPTypes);

numCTNodes=zeros(length(cPTypesU), numSIDT);
% nodeSID=G.Nodes.label(:, 2);

for ict=1:length(cPTypesU)
    if ~any(cPTypesU(ict)==fixedTypes)
    ctNodesi=cPTypes==cPTypesU(ict);
    ctNodesSIDi=nodeSID(ctNodesi);
    [ctNodesSIDiU, ~, jSIDU]=unique(ctNodesSIDi);
    cSIDi=accumarray(jSIDU, 1);
    UId=ismember(nodeSIDT,ctNodesSIDiU);

    numCTNodes(ict, UId)=cSIDi;
    end
end

numCTNodes=numCTNodes./sum(numCTNodes);