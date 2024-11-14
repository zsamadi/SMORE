function [nodeList, listLength, pathSections, subPathWC]=URPENH(gStruct, options)

vertexVec=(1:length(gStruct.neighs));
numVertex=length(vertexVec);
subPathCell=cell(numVertex,1);

gStructNeighs=gStruct.neighs;
listLength=[];

seedsAll=options.seedsAll;
iSelSeeds=~isempty(seedsAll);
ctypes=options.ctypes;
if iSelSeeds
    startEnd=seedsAll(:, [1, end]);
    startEnd=startEnd(:);
else
    startEnd=[];
end


coordinates=gStruct.coordinates;


parfor v=1:length(gStructNeighs)

    if iSelSeeds

        selectFlag=ismember(ctypes(v), startEnd);
    else
        selectFlag=true;
    end

    if selectFlag
        VExtension=gStructNeighs{v};
        subUPathMat=extendSubupath(v,VExtension,gStructNeighs,[], options);
        if iSelSeeds && ~isempty(subUPathMat)
            selectFlag=ismember(ctypes(subUPathMat), seedsAll, 'rows');
            subUPathMat=subUPathMat(selectFlag, :);
        end

    else
        subUPathMat=[];
    end

    subPathCell{v} =subUPathMat;

end
subPathVec=vertcat(subPathCell{:});
isRadialV=true(size(subPathVec, 1), 1);

if options.isChRadial

    info.coordinates=coordinates;
    for isp=1:size(subPathVec)

        isRadialV(isp)=checkRadial(subPathVec(isp, :), info);
    end
    subPathVec=subPathVec(isRadialV, :);
end

nodeList=num2cell(subPathVec, 2);
pathSections=ones(size(subPathVec, 1), 1);
subPathWC=ones(size(subPathVec, 1), 1);


function subUPathMat=extendSubupath(VSubgraph,VExtension, gStructNeighs,subUPathMat,options)
k=options.k;
pdv=options.pdv;

kSub=length(VSubgraph);
if kSub==k
    if VSubgraph(k)>VSubgraph(1)
        subUPathMat=[subUPathMat; VSubgraph];
    end
else

    while ~isempty(VExtension)
        z=VExtension(end);
        VExtension(end)=[];
        effectProb=pdv(length(VSubgraph)+1);
        randi=rand;
        goEhead=randi<effectProb;
        if (goEhead)
            zNeigh=gStructNeighs{z};

            if options.isChRadial
                VSubgraphNeigh=gStructNeighs(VSubgraph);
                VSubgraphNeigh=vertcat(VSubgraphNeigh{:});
                VSubgraphNeigh=unique(VSubgraphNeigh);
            else
                VSubgraphNeigh=[];
            end
            VExtensionN=setdiff(zNeigh, [VSubgraph(:);VSubgraphNeigh]);
            subUPathMat=extendSubupath([VSubgraph, z],VExtensionN, gStructNeighs, subUPathMat, options);
        end
    end

end
