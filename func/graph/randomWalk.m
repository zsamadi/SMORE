
function nodeList=randomWalk(gStruct, pathLength, numPaths)

alpha=0;


    numNodes = length(gStruct.neighs);

    nodeList = zeros(pathLength, numNodes, numPaths);

    for i =1:numPaths
        for j =1:numNodes
            indexList = generateSequence(j, gStruct, pathLength, alpha);

            nodeList(:,j, i)=indexList;

        end
    end
    nodeList=(reshape(nodeList,pathLength, [])).';
    isRadial=checkRadial(nodeList, gStruct);
    nodeList=nodeList(isRadial, :);

    
end



function result=generateSequence(startIndex, gStruct, pathLength, alpha)
    result=zeros(pathLength, 1);
    result(1) = startIndex;
    current = startIndex;

    randPathLength=rand(pathLength, 1);

    for i =2: pathLength
        if randPathLength(i) < alpha
            nextIndex = startIndex;
        else

            indexes=gStruct.neighs{current};
            probs = gStruct.weights{current};
            
            nextIndex =randsample(length(probs),1,true,probs) ;
            nextIndex=indexes(nextIndex);
            
        end
        result(i)= nextIndex;
        current = nextIndex;
    end
end






