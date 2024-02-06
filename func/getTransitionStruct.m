function gStruct=getTransitionStruct(network, isENodes)

if nargin<2
   isENodes=false;
end

if isENodes
    endNodes=network;
    weights=ones(size(endNodes, 1), 1);   
else
    endNodes=network.Edges.EndNodes;
    % weights=sqrt(network.Edges.Weight/min(network.Edges.Weight));
    % weights=network.Edges.Weight/median(network.Edges.Weight);
    weights=network.Edges.Weight;
    weights=weights/median(weights);

    weightsMin=min(weights);

    weightsMax=max(weights);
    weightsMd=median(weights);


end
numNodes=max(endNodes(:));

    neighsCell = cell(numNodes,1);
    weightsCell= cell(numNodes,1);




    for i =1: numNodes
        if isENodes
        neighsFlag= any(endNodes==i, 2);
        neighs =endNodes(neighsFlag, :);
        neighs=neighs.';
        neighs=neighs(:);
        neighs=neighs(neighs~=i);
        neighsCell{i}=neighs;
        weightsCell{i}=weights(neighsFlag);
        else
            neighsCell{i} = neighbors(network,i);
            eid = outedges(network,i);
            weightsCell{i}=weights(eid);
        end

    end

%     for i =1: numNodes
%         neighsWeights=weights(neighsFlagCell{i});
%         weightsCell{i}=neighsWeights/sum(neighsWeights);
%     end
    gStruct.neighs=neighsCell;
    gStruct.weights=weightsCell;
    
if isENodes
    gStruct.coordinates=zeros(length(gStruct.neighs), 2);
    gStruct.labels=ones(length(gStruct.neighs), 2);
 
else
    gStruct.coordinates=network.Nodes.Coordinates;
    gStruct.labels=network.Nodes.label;
end


gStruct.minWeight=weightsMin;
gStruct.maxWeight=weightsMax;
gStruct.mdWeight=weightsMd;









