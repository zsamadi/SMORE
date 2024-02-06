
function [G,gStruct,  h]=creatGraph(folderName, options)

numFileSection=options.retinaAndSection;
plotGraph=options.plotGraph;
isUniformWeight=options.isUniformWeight;

edgesUniqueTotal=[];
numcellsTotal=0;
xcoordsTotal=[];
ycoordsTotal=[];
cellSubtypeVecTotal=[];
WeigthsTotal=[];
cellSectionVecTotal=[];

if options.isRandom
      numCTypesRand=options.nCTypes;
       numNodesRand=options.nNodesRand;
       numNodesHold=floor(numNodesRand/5);
      
        sectNumberIdx=[1, 2];
        for sectNumberi=sectNumberIdx
            yshift=(sectNumberi-1)*sqrt(numNodesRand);

            if sectNumberi==1
                numNodesRandi=numNodesHold;
            else
                numNodesRandi=numNodesRand;
            end              
    
            xcoords=sqrt(numNodesRandi)*rand(numNodesRandi,1);
            xcoords=xcoords-mean(xcoords);
            ycoords=sqrt(numNodesRandi)*rand(numNodesRandi,1);
            ycoords=ycoords-mean(ycoords);
            
            cellSubtypeVec=randsample(numCTypesRand,numNodesRandi, true);
            cellSectionVec=ones(length(cellSubtypeVec),1);
         
            cellSectionVec(:)=sectNumberi;
    

    
            xcoordsTotal=[xcoordsTotal;xcoords];
            ycoordsTotal=[ycoordsTotal;ycoords+yshift];
            cellSubtypeVecTotal=[cellSubtypeVecTotal;cellSubtypeVec];
            cellSectionVecTotal=[cellSectionVecTotal;cellSectionVec];

    
            numcells=length(cellSubtypeVec);
    
    
            zcoords=[xcoords, ycoords];
    
%             Dmat = pdist2(zcoords,zcoords, 'euclidean');

            if strcmpi(options.gMode,"delaunay")
    
    
                Delauney_Triangle=delaunayTriangulation(zcoords);
      
                connectTriangle=Delauney_Triangle.ConnectivityList;
        
                edges=[connectTriangle(:,1:2); connectTriangle(:,2:3); connectTriangle(:,[1,3])];
            elseif strcmpi(options.gMode,"knn")

                numNeighs=options.numNeighs;

                nearNeighs = knnsearch(zcoords,zcoords, 'K',numNeighs);

                edges1=repelem(nearNeighs(:, 1), numNeighs-1, 1);
                edges2=nearNeighs(:, 2:end);
                edges2=edges2.';
                edges2=edges2(:);
                edges=[edges1, edges2];
            else
                % numNeighs=options.numNeighs;
                % [~, D] = knnsearch(zcoords,zcoords, 'K',numNeighs);
                rEps=options.rFEps;
                nearNeighs = rangesearch(zcoords,zcoords, rEps);
                edges=cell(length(nearNeighs),1);
                for iin=1:length(nearNeighs)
                    nearNeighsi=nearNeighs{iin};
                    numNeighs=length(nearNeighsi);
                    if numNeighs>1
                    edges1=repelem(iin, numNeighs-1, 1);
                    edges2=nearNeighsi(2:end);
                    edges2=edges2(:);
                    edges{iin}=[edges1, edges2];
                    else
                        edges{iin}=[];
                    end

                end
                edges=vertcat(edges{:});


            end



            edges=sort(edges, 2);
            [edgesUnique,~,~]=unique(edges, 'rows', 'stable');
    
%             edge_dist_idx=(edgesUnique(:,2)-1)*length(Dmat)+edgesUnique(:,1);
%             edge_dist=Dmat(edge_dist_idx);

            edge_dist02=zcoords(edgesUnique(:,2), :)-zcoords(edgesUnique(:,1), :);
            edge_dist02=edge_dist02.^2;
            edge_dist02=sum(edge_dist02,2);
            Weigths=exp(-edge_dist02/mean(sqrt(edge_dist02))^2);




    
            
          if strcmpi(options.gMode,"delaunay")
 
            weightThreshold=1e-6;
          else
              weightThreshold=0;
          end
          
    
            edgesUnique=edgesUnique(Weigths>weightThreshold, :);
            WeigthsValid=Weigths(Weigths>weightThreshold);
            edgesUniqueTotal=[edgesUniqueTotal;edgesUnique+numcellsTotal];
            numcellsTotal=numcellsTotal+numcells;
            WeigthsTotal=[WeigthsTotal;WeigthsValid];
        end
    cellSubtypeVec=cellSubtypeVecTotal;
    if isUniformWeight
        Weigths=ones(size(edgesUniqueTotal, 1),1);
    else
        Weigths=WeigthsTotal;
    end
    
    
    
    G=graph(edgesUniqueTotal(:,1),edgesUniqueTotal(:,2), Weigths);
    
    G.Nodes.label=[cellSubtypeVec, cellSectionVecTotal];
    G.Nodes.Coordinates=([xcoordsTotal, ycoordsTotal]);

% 
%             xcoords=rand(numNodesRand,1);
%             xcoords=numNodesRand*(xcoords-mean(xcoords));
%             ycoords=rand(numNodesRand,1);
%             ycoords=numNodesRand*(ycoords-mean(ycoords));
%             cellSubtypeVec=repmat((1:numCTypesRand).', numNodesRand/numCTypesRand,1);
%          
%             zcoords=[xcoords, ycoords];
%     
%             Dmat = pdist2(zcoords,zcoords, 'euclidean');
%     
%             % Delauney_Triangle = delaunay(xcoords,ycoords);
%     
%             Delauney_Triangle=delaunayTriangulation(zcoords);
%             connectTriangle=Delauney_Triangle.ConnectivityList;
%     
%             edges=[connectTriangle(:,1:2); connectTriangle(:,2:3); connectTriangle(:,[1,3])];
%             edges=sort(edges, 2);
%             [edgesUnique,~,~]=unique(edges, 'rows', 'stable');
%     
%             edge_dist_idx=(edgesUnique(:,2)-1)*length(Dmat)+edgesUnique(:,1);
%             % edge distances are extracted out from the computed distance matrix
%             edge_dist=Dmat(edge_dist_idx);
%             Weigths=exp(-edge_dist.^2/mean(edge_dist)^2);
%     
%             weightThreshold=1e-6;
%     
%             edgesUnique=edgesUnique(Weigths>weightThreshold, :);
%             WeigthsValid=Weigths(Weigths>weightThreshold);
% 
%     if isUniformWeight
%         Weigths=ones(size(WeigthsValid, 1),1);
%     else
%          Weigths=WeigthsValid;
%     end
%     
%     
%     G=graph(edgesUnique(:,1),edgesUnique(:,2), Weigths);
%     
%     G.Nodes.label=[cellSubtypeVec, ones(size(cellSubtypeVec))];
%     G.Nodes.Coordinates=zcoords;
%     xcoordsTotal=xcoords;
%     ycoordsTotal=ycoords;

else



    
    for ifile=1:numFileSection(1)
        csvName=strcat('Retina', num2str(ifile));
        filename=strcat(folderName,'\',csvName, '.csv');
        % filename='..\data\Retina1.csv';
    
        T = readtable(filename);
        sectNumber=T.SectionNumber;
    
        % cellSubtypeVec=cellSubtypeVec(randperm(length(cellSubtypeVec)));
    
        % iterMax=100;
        % s_cellSubtypeVec=zeros(numcells, iterMax);
        % for ii =1:iterMax
        %     s_cellSubtypeVec(:,ii)=cellSubtypeVec(randperm(numcells));
        % end

        yshift=(ifile-1)*10000;
    
    
        sectNumberIdx=(unique(sectNumber)).';
        sectNumberIdx=sectNumberIdx(~isnan(sectNumberIdx));
        secIDX=min(numFileSection(2),length(sectNumberIdx));
        for sectNumberi=sectNumberIdx(1:secIDX)
            xshift=(sectNumberi-1)*15000;
    
            xcoords=T.Remapped_X(sectNumber==sectNumberi);
            xcoords=xcoords-mean(xcoords);

            if options.randCoords
                xcoordsMin=min(xcoords);
                xcoordsMax=max(xcoords);
                xcoords=(xcoordsMax-xcoordsMin)*rand(length(xcoords),1)+xcoordsMin;
            end

            



            ycoords=T.Remapped_Y(sectNumber==sectNumberi);
            ycoords=ycoords-mean(ycoords);

            if options.randCoords          
                ycoordsMin=min(ycoords);
                ycoordsMax=max(ycoords);
                ycoords=(ycoordsMax-ycoordsMin)*rand(length(ycoords),1)+ycoordsMin;
            end

            cellSubtypeVec=T.Subtype(sectNumber==sectNumberi);
            cellSectionVec=ones(length(cellSubtypeVec),1);
         
            cellSectionVec(:)=(ifile-1)*secIDX+sectNumberi;
    
            %
            %         xcoords=xcoords(1:5);
            %         ycoords=ycoords(1:5);
            %         cellSubtypeVec=cellSubtypeVec(1:5);
    
            xcoordsTotal=[xcoordsTotal;xcoords+xshift];
            ycoordsTotal=[ycoordsTotal;ycoords+yshift];
            cellSubtypeVecTotal=[cellSubtypeVecTotal;cellSubtypeVec];
            cellSectionVecTotal=[cellSectionVecTotal;cellSectionVec];
    
            % xcoords=xcoords(cellSubtypeVec~=15);
            % ycoords=ycoords(cellSubtypeVec~=15);
            % cellSubtypeVec=cellSubtypeVec(cellSubtypeVec~=15);
    
            % %
    
            numcells=length(cellSubtypeVec);
    
    
            zcoords=[xcoords, ycoords];
    
%             Dmat = pdist2(zcoords,zcoords, 'euclidean');
    
            % Delauney_Triangle = delaunay(xcoords,ycoords);
    


            if strcmpi(options.gMode,"delaunay")
    
    
                Delauney_Triangle=delaunayTriangulation(zcoords);
      
                connectTriangle=Delauney_Triangle.ConnectivityList;
        
                edges=[connectTriangle(:,1:2); connectTriangle(:,2:3); connectTriangle(:,[1,3])];
            elseif strcmpi(options.gMode,"knn")

                numNeighs=options.numNeighs;

                nearNeighs = knnsearch(zcoords,zcoords, 'K',numNeighs);

                edges1=repelem(nearNeighs(:, 1), numNeighs-1, 1);
                edges2=nearNeighs(:, 2:end);
                edges2=edges2.';
                edges2=edges2(:);
                edges=[edges1, edges2];
            else
                numNeighs=2;
                [~, D] = knnsearch(zcoords,zcoords, 'K',numNeighs);
                DS=sort(D(:));
                rEps=DS(end-10);
                nearNeighs = rangesearch(zcoords,zcoords, rEps);
                edges=cell(length(nearNeighs),1);
                for iin=1:length(nearNeighs)
                    nearNeighsi=nearNeighs{iin};
                    numNeighs=length(nearNeighsi);
                    if numNeighs>1
                    edges1=repelem(iin, numNeighs-1, 1);
                    edges2=nearNeighsi(2:end);
                    edges2=edges2(:);
                    edges{iin}=[edges1, edges2];
                    else
                        edges{iin}=[];
                    end

                end
                edges=vertcat(edges{:});


            end
    
            % figure
            % triplot(Delauney_Triangle)

    
            edges=sort(edges, 2);
            [edgesUnique,~,~]=unique(edges, 'rows', 'stable');
    
%             edge_dist_idx=(edgesUnique(:,2)-1)*length(Dmat)+edgesUnique(:,1);
            % edge distances are extracted out from the computed distance matrix
%             edge_dist=Dmat(edge_dist_idx);
%             Weigths=exp(-edge_dist.^2/mean(edge_dist)^2);


            edge_dist02=zcoords(edgesUnique(:,2), :)-zcoords(edgesUnique(:,1), :);
            edge_dist02=edge_dist02.^2;
            edge_dist02=sum(edge_dist02,2);
            Weigths=exp(-edge_dist02/mean(sqrt(edge_dist02))^2);

          if strcmpi(options.gMode,"delaunay")
 
            weightThreshold=1e-6;
          else
              weightThreshold=0;
          end

    
            edgesUnique=edgesUnique(Weigths>weightThreshold, :);
            WeigthsValid=Weigths(Weigths>weightThreshold);
            edgesUniqueTotal=[edgesUniqueTotal;edgesUnique+numcellsTotal];
            numcellsTotal=numcellsTotal+numcells;
            WeigthsTotal=[WeigthsTotal;WeigthsValid];
        end
    end
    cellSubtypeVec=cellSubtypeVecTotal;
    if isUniformWeight
        Weigths=ones(size(edgesUniqueTotal, 1),1);
    else
        Weigths=WeigthsTotal;
    end
    
    
    
    G=graph(edgesUniqueTotal(:,1),edgesUniqueTotal(:,2), Weigths);
    
    G.Nodes.label=[cellSubtypeVec, cellSectionVecTotal];
    G.Nodes.Coordinates=([xcoordsTotal, ycoordsTotal]);

end

gStruct=getTransitionStruct(G);
% gStruct.coordinates=G.Nodes.Coordinates;
% gStruct.labels=G.Nodes.label;

if plotGraph
    figure
    h=plot(G,'XData',xcoordsTotal,'YData',ycoordsTotal);
    grid on
    hold on
else
    h=0;
end