function [nodeList, listLength, pathSections, subPathWC]=enumerateUSRKPathsTest(gStruct, options)
%created on 230902 to test if checking for radiality in mid subgraf lengths helps to improve speed 

 k=options.k;
 pdv=options.pdv;



% gStruct.neighs={4, 4, 4, [1, 2, 3, 12, 8], 8, 8, 8, [5, 6, 7, 12, 4], 12, 12, 12, [9, 10, 11, 4, 8]};
% gStruct=gStruct.neighs;

% pdv=pdv(2:end)/gStruct.minWeight;
pdv=pdv(2:end);

% pdvMax=1/gStruct.minWeight;
% gStruct.pdvMax=pdvMax;

% pdv=pdv/gStruct.minWeight;
% 
% 

% pdv=[pdv(:); pdvM];

vertexVec=(1:length(gStruct.neighs));
numVertex=length(vertexVec);
subPathCell=cell(numVertex,1);
subPathWC=cell(numVertex,1);
pathSections=cell(numVertex,1);
gStructNeighs=gStruct.neighs;

gStructWeights=gStruct.weights;
options.pdv=pdv;


% sampVertex=(1:numVertex);
% 
% sampVertex(notSampleVertex)=[];
listLength=zeros(length(gStructNeighs), 1);

parfor v=1:length(gStructNeighs)
    % v=sampVertex(iv);



    VExtension=gStructNeighs{v};
    WExtension=gStructWeights{v};

    % VExtension=repelem(VExtension, gStructWeights{v});

    [subUPathMat, wSubG]=extendSubupath(v,VExtension,WExtension,gStruct,[],1, [], options);
    % tic

    subPathCell{v} =subUPathMat;
    subPathWC{v} =wSubG;

    pathSections{v}=repelem(gStruct.labels(v, 2), size(subPathCell{v}, 1),1);
    listLength(v)=size(subUPathMat, 1);
    % toc

    % if (mod(v, 500)==1)
    %     fprintf('%d vertices processed so far \n', v);
    % end
end
subPathVec=vertcat(subPathCell{:});
subPathWC=vertcat(subPathWC{:});

nodeList=num2cell(subPathVec, 2);

pathSections=vertcat(pathSections{:});

% listLength=repelem(k, length(pathSections), 1);



function [subUPathMat, wSubG]=extendSubupath(VSubgraph,VExtension,WExtension, gStruct,subUPathMat, wGoAhead, wSubG, options)
 k=options.k;
 pdv=options.pdv;

kSub=length(VSubgraph);
if kSub==k
  if VSubgraph(k)>VSubgraph(1)
%      rvDist=gStruct.coordinates(VSubgraph(1:k-1), :)-gStruct.coordinates(VSubgraph(k), :);
%      rvDist=sum(rvDist.^2, 2);

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

%      if ~goEhead && kSub>7
%          if VSubgraph(end)>v
%             subUPathCell=[subUPathCell,0, VSubgraph];
%          end
%              
%      end

    while ~isempty(VExtension)

        z=VExtension(end);
        % zFalgs=VExtension==z;
        % numChances=sum(zFalgs);
        VExtension(end)=[];
        w=WExtension(end);
        WExtension(end)=[];

        effectProb=pdv(length(VSubgraph));

            % randChance=rand(1, numChances);
            % goAheads=randChance<pdv(length(VSubgraph));

            randi=rand;

            goEhead=randi<effectProb;
        
        
        
        if (goEhead)
            % zdistn=norm(gStruct.coordinates(z, :)-gStruct.coordinates(v, :));
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
                % VExtensionN=repelem(VExtensionN, zWeights(iz));

                [subUPathMat, wSubG]=extendSubupath([VSubgraph, z],VExtensionN,WExtensionN , gStruct, subUPathMat, wGoAheadN, wSubG, options);
        end
    end

end
