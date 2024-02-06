
function [G,gStruct,  h, BregmaAllU, AnimalIDIdx, TOut, edgeDistStd]=creatHPGraph(T,folderName, options)

numIDs=options.hypoIDs;
plotGraph=options.plotGraph;
isUniformWeight=options.isUniformWeight;

edgesUniqueTotal=[];
numcellsTotal=0;
xcoordsTotal=[];
ycoordsTotal=[];
cellSubtypeVecTotal=[];
WeigthsTotal=[];
cellSectionVecTotal=[];




% filename='..\data\Retina1.csv';


AnimalID=T.Animal_ID;

% cellSubtypeVec=cellSubtypeVec(randperm(length(cellSubtypeVec)));

% iterMax=100;
% s_cellSubtypeVec=zeros(numcells, iterMax);
% for ii =1:iterMax
%     s_cellSubtypeVec(:,ii)=cellSubtypeVec(randperm(numcells));
% end


BregmaAll=T.Bregma;
BregmaAllU=unique(BregmaAll);
% BregmaAllIndex=(1:length(BregmaAllU));


AnimalIDIdx=(unique(AnimalID)).';
AnimalIDIdx=AnimalIDIdx(~isnan(AnimalIDIdx));
secIDX=min(numIDs,length(AnimalIDIdx));
AnimalIDIdx=AnimalIDIdx(1:secIDX);
% figure
TCia=cell(length(AnimalIDIdx), 1);
edgeDistSum=[0,0];
if options.xyNoiseStd>0
    cellTAll=T.cellSubtype;
    notFixedTypes=~ismember(cellTAll, options.fixedTypes);

    if options.isNoiseCTS



        cellTAll(notFixedTypes)=cellTAll(notFixedTypes)+rand(sum(notFixedTypes),1);
        [~, iu, ju]=unique(cellTAll);
        xyNoiseU=options.xyNoiseStd*(rand(length(iu), 2)-0.5);
        xyNoise=xyNoiseU(ju, :);


    else

        xyNoise=zeros(size(T, 1), 2);
        xyNoiseNF=rand(sum(notFixedTypes), 2);
        xyNoiseNF=(xyNoiseNF-mean(xyNoiseNF))./std(xyNoiseNF);
        xyNoiseNF=options.xyNoiseStd*xyNoiseNF;
        xyNoise(notFixedTypes, :)=xyNoiseNF;



    end
else
    xyNoise=zeros(size(T, 1), 2);
end
numNodes=0;

nearNeighsia=cell(length(AnimalIDIdx),1);

for iAnimalID=1:length(AnimalIDIdx)
    AnimalIDi=AnimalIDIdx(iAnimalID);
    loopFlafID=(AnimalID==AnimalIDi);
    T1=T(loopFlafID, :);
    xyNoise1=xyNoise(loopFlafID, :);

    Bregma=T1.Bregma;
    BregmaU=unique(Bregma);
    if length(BregmaU)<2
        yshift=(iAnimalID-1)*3000;
    else
        yshift=(AnimalIDi-1)*3000;
    end
    [~, BregmasInAll]=ismember(BregmaU, BregmaAllU);
    TCib=cell(length(BregmaU), 1);
    nearNeighsib=cell(length(BregmaU), 1);

    % for iBregma=1:min(1,length(BregmaU))
    for iBregma=1:length(BregmaU)

       

        xshift=(BregmaU(iBregma))*1e5/2;
        loopFlagBrg=(Bregma==BregmaU(iBregma));
        T2=T1(loopFlagBrg, :);
        xyNoise2=xyNoise1(loopFlagBrg, :);

        % randNoiseIdx=randi(size(xyNoise, 1), size(T2, 1), 1);

        maxX=max(T2.Centroid_X);
        minX=min(T2.Centroid_X);
       

        xcoords=T2.Centroid_X+xyNoise2(:, 1);



        gFlag=xcoords>maxX;
        lFlag=xcoords<minX;
       

        xcoords(gFlag|lFlag)=xcoords(gFlag|lFlag)-0.9*xyNoise2(gFlag|lFlag, 1);
        % xcoords(lFlag)=xcoords(lFlag)-xyNoise2(lFlag, 1);

        xcoords=xcoords-mean(xcoords);


        maxY=max(T2.Centroid_Y);
        minY=min(T2.Centroid_Y);


        ycoords=T2.Centroid_Y+xyNoise2(:, 2);

        gFlag=ycoords>maxY;
        lFlag=ycoords<minY;

        ycoords(gFlag|lFlag)=ycoords(gFlag|lFlag)-0.9*xyNoise2(gFlag|lFlag, 2);
        % ycoords(lFlag)=ycoords(lFlag)-xyNoise2(lFlag, 2);


        ycoords=ycoords-mean(ycoords);
        % xyNoise(1:sum(loopFlagBrg), :)=[];

        % Neuron_cluster_ID


        cellSubtypeVec=T2.cellSubtype;

        cellSectionVec=ones(length(cellSubtypeVec),1);

        cellSectionVec(:)=(AnimalIDi-1)*length(BregmaAllU)+BregmasInAll(iBregma);

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
        

        cellSubtypeVecU=unique(cellSubtypeVec);
        % figure
        % hold on
        % for ii=1:length(cellSubtypeVecU)
        %     zcoordsi=zcoords(cellSubtypeVec==cellSubtypeVecU(ii), :);
        %     scatter(zcoordsi(:, 1), zcoordsi(:, 2), 9, 'filled')
        % end
        % plot([200; 700], [-800; -800], '-k', 'LineWidth',2)
        % 
        % text(450,-720, '\bf 500 \mum', 'HorizontalAlignment','center','FontSize',16)
        check=1;











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

           rEps=options.rEpsilon;
           neighsi= rangesearch(zcoords,zcoords, rEps);
    
           if numNodes>0
                for iNd=1:length(neighsi)
                     neighsi{iNd}=neighsi{iNd}+numNodes;
                end
           end
           numNodesi=length(xcoords);
           numNodes=numNodes+numNodesi;
           
           nearNeighsib{iBregma}=neighsi;
       

        % figure
        % triplot(Delauney_Triangle)


        edges=sort(edges, 2);
        [edgesUnique,~,~]=unique(edges, 'rows', 'stable');

        %             edge_dist_idx=(edgesUnique(:,2)-1)*length(Dmat)+edgesUnique(:,1);
        % edge distances are extracted out from the computed distance matrix
        %             edge_dist=Dmat(edge_dist_idx);
        %             Weigths=exp(-edge_dist.^2/mean(edge_dist)^2);


        edgeDist02=zcoords(edgesUnique(:,2), :)-zcoords(edgesUnique(:,1), :);
        edgeDist02=edgeDist02.^2;
        edgeDist02=sum(edgeDist02,2);
        Weigths=exp(-edgeDist02/mean(sqrt(edgeDist02))^2);



        

        weightThreshold=1e-10;
        edgeDistSum=edgeDistSum+[sum(sqrt(edgeDist02(Weigths>weightThreshold))), sum(Weigths>weightThreshold)];



        edgesUnique=edgesUnique(Weigths>weightThreshold, :);
        WeigthsValid=Weigths(Weigths>weightThreshold);
        edgesUniqueTotal=[edgesUniqueTotal;edgesUnique+numcellsTotal];
        numcellsTotal=numcellsTotal+numcells;
        WeigthsTotal=[WeigthsTotal;WeigthsValid];

        % 
        % G=graph(edgesUniqueTotal(:,1),edgesUniqueTotal(:,2), WeigthsTotal);
        % 
        % G.Nodes.Coordinates=([xcoordsTotal, ycoordsTotal]);
        % plot(G,'XData',xcoordsTotal,'YData',ycoordsTotal);
        % grid on
        % hold on
        TCib{iBregma}=T2;




    end

    TCia{iAnimalID} = vertcat(TCib{:});
    nearNeighsia{iAnimalID}=vertcat(nearNeighsib{:});
end

TOut=vertcat(TCia{:});

TOut.Centroid_X=xcoordsTotal;
TOut.Centroid_Y=ycoordsTotal;


cellSubtypeVec=cellSubtypeVecTotal;
if isUniformWeight
    Weigths=ones(size(edgesUniqueTotal, 1),1);
else
    Weigths=WeigthsTotal;
end

nearNeighs=vertcat(nearNeighsia{:});

G=graph(edgesUniqueTotal(:,1),edgesUniqueTotal(:,2), Weigths);

G.Nodes.label=[cellSubtypeVec, cellSectionVecTotal];
G.Nodes.Coordinates=([xcoordsTotal, ycoordsTotal]);
if ~isempty(nearNeighs)
    G.Nodes.nearNeighs=nearNeighs;
end

[cellSubtypeVecU, ~, jU]=unique(cellSubtypeVec);

cellCnts=accumarray(jU, 1);
[~, ist]=sort(cellCnts, 'descend');

cellSubtypeVecU=cellSubtypeVecU(ist);


gStruct=getTransitionStruct(G);
% gStruct.coordinates=G.Nodes.Coordinates;
% gStruct.labels=G.Nodes.label;


nodesAll=(1:length(cellSubtypeVec));

edgeDistStd=edgeDistSum(1)/edgeDistSum(2);

if plotGraph

    rng(1)

    numHighlight=min(10, length(cellSubtypeVecU));
    % figure('WindowState' ,'maximize')
    figure
    ncolor=jet(256);
    % ncolor=ncolor(end:-1:1, :);
    ncolor=ncolor(randperm(256), :);
    ncolor=ncolor(1:numHighlight, :);

    h=plot(G,'XData',xcoordsTotal,'YData',ycoordsTotal, 'EdgeColor',[0,0,0], 'MarkerSize',1);
    for iu=1:length(cellSubtypeVecU)
        % ncolor=rand(1, 3);

        if iu<=numHighlight
            highlight(h,nodesAll(cellSubtypeVec==cellSubtypeVecU(iu)),'NodeColor',ncolor(iu, :), 'EdgeColor',[1,1,1])
        else
            highlight(h,nodesAll(cellSubtypeVec==cellSubtypeVecU(iu)),'NodeColor',[0.6980    0.7451    0.7098], 'EdgeColor',[1,1,1])
        end


    end
        hold on
        ax=zeros(numHighlight+1, 1);

    for iu=1:numHighlight
        ax(iu)=plot(NaN,NaN,'.', 'MarkerFaceColor', ncolor(iu, :), 'MarkerEdgeColor', ncolor(iu, :), 'markersize', 20); %plotting invisible points of desired colors
    end 
        ax(end)=plot(NaN,NaN,'.', 'MarkerFaceColor', [0.6980    0.7451    0.7098], 'MarkerEdgeColor',[0.6980    0.7451    0.7098], 'markersize', 20); %plotting invisible points of desired colors


    legendText=lower(options.cellTypesOne(cellSubtypeVecU(1:numHighlight)));
    legendText=[legendText;"others"];

    % legendText=num2str(cellSubtypeVecU(1:numHighlight));

    % legend(ax,legendText,'Orientation','horizontal')
    legend(ax,legendText,'NumColumns',1, 'Location', 'Bestoutside')


    set(gca,'xtick',BregmaAllU*1e5/2,'xticklabel',BregmaAllU) 
    if length(BregmaAllU)<2
        set(gca,'ytick',(0:length(AnimalIDIdx)-1)*3000,'yticklabel',AnimalIDIdx)
    else
        set(gca,'ytick',(AnimalIDIdx-1)*3000,'yticklabel',AnimalIDIdx)
    end
    xlabel('bregma')
    ylabel('animal\_ID')
    
    grid minor

    if options.xyNoiseStd>0
         figname=strcat(folderName, 'NoisyGraphWithTypes.png');

    else

        figname=strcat(folderName, 'GraphWithTypes.png');
    end
    % print(gcf,figname(1:end-4), '-djpeg', '-r600'); %<-Save as jpg with 600 DPI
    % print(gcf,figname(1:end-4), '-djpeg'); %<-Save as jpg 

    saveas(gcf,figname)





    % cellTypesOne=vertcat(cellTypesOut(1:14, 2), cellTypesOut(15:end, 1));
    % cPTypes=(G.Nodes.label(:, 1)).';
    % cPTypesU=unique(cPTypes);

    % PWMT=sum(cPTypes.'==cPTypesU);
    % PWMT=PWMT/sum(PWMT);
    % numFixedTypes=length(PWMT);
    % [~, fixedTypes]=sort(PWMT, 'descend');
    %     fixedTypes=fixedTypes(1:numFixedTypes);
    %     fixedTypes=fixedTypes(:);
        % 
        % for iFx=1:length(fixedTypes)
        % 
        %     fig=figure('WindowState' ,'maximize');
        % 
        %     h=plot(G,'XData',xcoordsTotal,'YData',ycoordsTotal);
        %         % ncolor=rand(1, 3);
        %     highlight(h,nodesAll(cellSubtypeVec==fixedTypes(iFx)),'NodeColor','g')
        %     hold on
        % 
        % 
        %    ax=plot(NaN,NaN,'.', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'markersize', 20); %plotting invisible points of desired colors
        % 
        % 
        %     legendText=cellTypesOne(fixedTypes(iFx));
        % 
        %     % legend(ax,legendText,'Orientation','horizontal')
        %     legend(ax,legendText)
        % 
        % 
        %     set(gca,'xtick',BregmaAllU*1e5/2,'xticklabel',BregmaAllU) 
        %     set(gca,'ytick',(AnimalIDIdx(1:secIDX)-1)*3000,'yticklabel',(1:secIDX))
        %     xlabel('Bregma')
        %     ylabel('Animal\_ID')
        % 
        %     grid minor
        %     figname=strcat('output\GraphWithSubtTypes' , num2str(fixedTypes(iFx)), '.fig');
        %     % print(gcf,figname(1:end-4), '-djpeg', '-r600'); %<-Save as jpg with 600 DPI
        %     print(gcf,figname(1:end-4), '-djpeg'); %<-Save as jpg 
        %     close(fig)
        % end




else
    h=0;
end

check=1;
