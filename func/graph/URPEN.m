function [nodeList, listLength, pathSections, subPathWC]=URPEN(gStruct, options)

vertexVec=(1:length(gStruct.neighs));
numVertex=length(vertexVec);
subPathCell=cell(numVertex,1);
subPathWC=cell(numVertex,1);
pathSections=cell(numVertex,1);
gStructNeighs=gStruct.neighs;
gStructWeights=gStruct.weights;
listLength=zeros(length(gStructNeighs), 1);

parfor v=1:length(gStructNeighs)

    VExtension=gStructNeighs{v};
    WExtension=gStructWeights{v};
    [subUPathMat, wSubG]=extendSubupath(v,VExtension,WExtension,gStruct,[],1, [], options);
    subPathCell{v} =subUPathMat;
    subPathWC{v} =wSubG;
    pathSections{v}=repelem(gStruct.labels(v, 2), size(subPathCell{v}, 1),1);
    listLength(v)=size(subUPathMat, 1);

end
subPathVec=vertcat(subPathCell{:});
subPathWC=vertcat(subPathWC{:});

nodeList=num2cell(subPathVec, 2);

pathSections=vertcat(pathSections{:});



function [subUPathMat, wSubG]=extendSubupath(VSubgraph,VExtension,WExtension, gStruct,subUPathMat, wGoAhead, wSubG, options)
k=options.k;
pdv=options.pdv;

kSub=length(VSubgraph);
if kSub==k
    if VSubgraph(k)>VSubgraph(1)
        if options.isChRadial

            if (checkRadial(VSubgraph, gStruct))
                subUPathMat=[subUPathMat; VSubgraph];
                wSubG=[wSubG;wGoAhead];

            end
        else
            subUPathMat=[subUPathMat; VSubgraph];
            wSubG=[wSubG;wGoAhead];
        end

    end
else

    while ~isempty(VExtension)
        z=VExtension(end);
        VExtension(end)=[];
        w=WExtension(end);
        WExtension(end)=[];
        effectProb=pdv(length(VSubgraph)+1);
        randi=rand;
        goEhead=randi<effectProb;
        if (goEhead)
            wGoAheadN=w*wGoAhead;
            zNeigh=gStruct.neighs{z};
            zWeights=gStruct.weights{z};
            if options.isChRadial
                VSubgraphNeigh=gStruct.neighs(VSubgraph);
                VSubgraphNeigh=vertcat(VSubgraphNeigh{:});
                VSubgraphNeigh=unique(VSubgraphNeigh);
            else
                VSubgraphNeigh=[];
            end

            [VExtensionN, iz]=setdiff(zNeigh, [VSubgraph(:);VSubgraphNeigh]);
            WExtensionN=zWeights(iz);
            [subUPathMat, wSubG]=extendSubupath([VSubgraph, z],VExtensionN,WExtensionN , gStruct, subUPathMat, wGoAheadN, wSubG, options);
        end
    end

end
