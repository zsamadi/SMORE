
function [G,gStruct,  cellSubtypeVecU]=creat3DGraph(csvName, folderName,options)

plotGraph=options.plotGraph;
isUniformWeight=options.isUniformWeight;


cellSectionVecTotal=[];

    
        filename=strcat(csvName, '.csv');
        % filename='..\data\Retina1.csv';
    
        T = readtable(filename);
    
        % cellSubtypeVec=cellSubtypeVec(randperm(length(cellSubtypeVec)));
    
        % iterMax=100;
        % s_cellSubtypeVec=zeros(numcells, iterMax);
        % for ii =1:iterMax
        %     s_cellSubtypeVec(:,ii)=cellSubtypeVec(randperm(numcells));
        % end
    
        xshift=0;
    
            yshift=0;

            isDoublet=T.doublet;

            isDoubletTF=false(length(isDoublet), 1);

            for iDBL=1:length(isDoublet)
                isDoubletTF(iDBL)=isDoublet{iDBL}=="TRUE";
            end

            T=T(~isDoubletTF, :);


    
            xcoords=T.global_x;
            xcoords=xcoords-mean(xcoords);

            ycoords=T.global_y;
            ycoords=ycoords-mean(ycoords);

            zcoords=T.global_z;
            zcoords=zcoords-mean(zcoords);

            cellSubtypeVec=T.cluster;
         
    
            %
            %         xcoords=xcoords(1:5);
            %         ycoords=ycoords(1:5);
            %         cellSubtypeVec=cellSubtypeVec(1:5);
    
            xcoordsTotal=xcoords+xshift;
            ycoordsTotal=ycoords+yshift;
            zcoordsTotal=zcoords;

    
            % xcoords=xcoords(cellSubtypeVec~=15);
            % ycoords=ycoords(cellSubtypeVec~=15);
            % cellSubtypeVec=cellSubtypeVec(cellSubtypeVec~=15);
    
            % %
    
            numcells=length(cellSubtypeVec);
    
    
            xyzcoords=[xcoords, ycoords, zcoords];
    
%             Dmat = pdist2(zcoords,zcoords, 'euclidean');
    
            % Delauney_Triangle = delaunay(xcoords,ycoords);
    


            if strcmpi(options.gMode,"delaunay")
    
    
                Delauney_Triangle=delaunayTriangulation(xyzcoords);
      
                connectTriangle=Delauney_Triangle.ConnectivityList;
        
                edges=[connectTriangle(:,[1,2]); connectTriangle(:,[1,3]); connectTriangle(:,[1,4]); connectTriangle(:,[2,3]); connectTriangle(:,[2,4]); connectTriangle(:,[3,4])];
            elseif strcmpi(options.gMode,"knn")

                numNeighs=options.numNeighs;

                nearNeighs = knnsearch(xyzcoords,xyzcoords, 'K',numNeighs);

                edges1=repelem(nearNeighs(:, 1), numNeighs-1, 1);
                edges2=nearNeighs(:, 2:end);
                edges2=edges2.';
                edges2=edges2(:);
                edges=[edges1, edges2];
            else
                numNeighs=2;
                [~, D] = knnsearch(xyzcoords,xyzcoords, 'K',numNeighs);
                DS=sort(D(:));
                rEps=DS(end-10);
                nearNeighs = rangesearch(xyzcoords,xyzcoords, rEps);
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


            edge_dist02=xyzcoords(edgesUnique(:,2), :)-xyzcoords(edgesUnique(:,1), :);
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
            edgesUniqueTotal=edgesUnique;
            WeigthsTotal=WeigthsValid;
  
    if isUniformWeight
        Weigths=ones(size(edgesUniqueTotal, 1),1);
    else
        Weigths=WeigthsTotal;
    end
    
    
    
    G=graph(edgesUniqueTotal(:,1),edgesUniqueTotal(:,2), Weigths);

    cell_type=T.cell_type;
    subclass=T.subclass;

    [cellSubtypeVecU, iU, jU]=unique(cellSubtypeVec);

    cell_typeU=cell_type(iU);
    nnIdx=zeros(length(cellSubtypeVecU), 1);
    iCSVn=1;
    for iCSV=1:length(cellSubtypeVecU)
        if strcmpi(cellSubtypeVecU{iCSV},cell_typeU{iCSV})
            nnIdx(iCSVn)=iCSV;
            iCSVn=iCSVn+1;
        else
            cellSubtypeVecU{iCSV}=strcat(cell_typeU{iCSV}(1),'-', cellSubtypeVecU{iCSV});
        end
    end
    nnIdx=nnIdx(nnIdx>0);
    sortIdx=(1:length(cellSubtypeVecU));
    sortIdx(nnIdx)=[];
    sortIdx=[nnIdx;sortIdx(:)];
    cellSubtypeVecU=cellSubtypeVecU(sortIdx);

    sortIdxInv=sortIdx;
    sortIdxInv(sortIdx)=(1:length(cellSubtypeVecU));


    jU=sortIdxInv(jU);



    rEps=options.rEpsilon;
    G.Nodes.nearNeighs= rangesearch(xyzcoords,xyzcoords, rEps);
    



    subCNTs=accumarray(jU, 1);
    [subCNTs, iSrt]=sort(subCNTs, 'descend');

    cellSubtypeVecd=(1:length(iU));
    cellSubtypeVecUDS=cellSubtypeVecd(iSrt);
    cellSubtypeVecUS=cellSubtypeVecU(iSrt);

    cellSubtypeVecd=cellSubtypeVecd(jU);

    G.Nodes.label=[cellSubtypeVecd(:), ones(numcells, 1)];
    G.Nodes.Coordinates=([xcoordsTotal, ycoordsTotal, zcoordsTotal]);


gStruct=getTransitionStruct(G);
% gStruct.coordinates=G.Nodes.Coordinates;
% gStruct.labels=G.Nodes.label;

if plotGraph

    nodesAll=(1:length(cellSubtypeVec));

    numHighlight=40;

        ncolor=jet(numHighlight+1);
    % ncolor=ncolor(end:-1:1, :);
    ncolor=ncolor(randperm(numHighlight+1), :);
    ncolor=ncolor(1:numHighlight, :);

    figure('WindowState','maximized')
    h=plot(G,'XData',xcoordsTotal,'YData',ycoordsTotal,'ZData',zcoordsTotal,'NodeColor',[0.8    0.8    0.8], 'EdgeColor',[0.8    0.8    0.8]);

    for iu=1:length(cellSubtypeVecUDS)
        % ncolor=rand(1, 3);

        if iu<=numHighlight
            highlight(h,nodesAll(cellSubtypeVecd==cellSubtypeVecUDS(iu)),'NodeColor',ncolor(iu, :), 'EdgeColor',[0.6    0.7451    0.7098])
        else
            highlight(h,nodesAll(cellSubtypeVecd==cellSubtypeVecUDS(iu)),'NodeColor',[0.6980    0.7451    0.7098], 'EdgeColor',[0.6980    0.7451    0.7098])
        end


    end
        hold on
        ax=zeros(numHighlight+1, 1);

    for iu=1:numHighlight
        ax(iu)=plot(NaN,NaN,'.', 'MarkerFaceColor', ncolor(iu, :), 'MarkerEdgeColor', ncolor(iu, :), 'markersize', 20); %plotting invisible points of desired colors
    end 
        ax(end)=plot(NaN,NaN,'.', 'MarkerFaceColor', [0.6980    0.7451    0.7098], 'MarkerEdgeColor',[0.6980    0.7451    0.7098], 'markersize', 20); %plotting invisible points of desired colors


    legendText=cellSubtypeVecUS(1:numHighlight);
    legendText=[legendText;"others"];

    % legendText=num2str(cellSubtypeVecU(1:numHighlight));

    % legend(ax,legendText,'Orientation','horizontal')
    legend(ax,legendText,'NumColumns',2, 'Location', 'Bestoutside')

    grid on
    hold on

    if options.xyNoiseStd>0
         figname=strcat(folderName, 'NoisyGraphWithTypes.png');
         % fignameMax=strcat(folderName, 'NoisyGraphWithTypesMax.png');
    else

        figname=strcat(folderName, 'GraphWithTypes.png');
        % fignameMax=strcat(folderName, 'GraphWithTypesMax.png');
   end


   saveas(gcf,figname)


else
    h=0;
end