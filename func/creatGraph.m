
function [G,gStruct,  cellTypes]=creatGraph(folderName, options)

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


else



    
    for ifile=1:numFileSection(1)
        csvName=strcat('Retina', num2str(ifile));
        filename=strcat(folderName,'\',csvName, '.csv');
    
        T = readtable(filename);
        sectNumber=T.SectionNumber;


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
    
    
            xcoordsTotal=[xcoordsTotal;xcoords+xshift];
            ycoordsTotal=[ycoordsTotal;ycoords+yshift];
            cellSubtypeVecTotal=[cellSubtypeVecTotal;cellSubtypeVec];
            cellSectionVecTotal=[cellSectionVecTotal;cellSectionVec];
    

            numcells=length(cellSubtypeVec);
    
    
            zcoords=[xcoords, ycoords];

    


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
    
    
            edges=sort(edges, 2);
            [edgesUnique,~,~]=unique(edges, 'rows', 'stable');
    


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
    cellTypes=unique(cellSubtypeVec);
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

if plotGraph
    figure
    h=plot(G,'XData',xcoordsTotal,'YData',ycoordsTotal);
    grid on
    hold on
else
    h=0;
end