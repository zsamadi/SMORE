
function [G,gStruct,  h, sectionU, AnimalIDIdx, TOut, edgeDistStd]=creatBipolarGraph(DTable,folderName, options)




plotGraph=options.plotGraph;
isUniformWeight=options.isUniformWeight;

edgesUniqueTotal=[];
numcellsTotal=0;
xcoordsTotal=[];
ycoordsTotal=[];
subtypeVecTotal=[];
WeigthsTotal=[];
cellSectionVecTotal=[];




% filename='..\data\Retina1.csv';


AnimalID=DTable.animalID;

% subtype_codesVec=subtype_codesVec(randperm(length(subtype_codesVec)));

% iterMax=100;
% s_subtype_codesVec=zeros(numcells, iterMax);
% for ii =1:iterMax
%     s_subtype_codesVec(:,ii)=subtype_codesVec(randperm(numcells));
% end


sectionAll=DTable.SectionNumber;
sectionU=unique(sectionAll);
% sectionid_codesAllIndex=(1:length(sectionid_codesAllU));


AnimalIDIdx=(unique(AnimalID)).';
AnimalIDIdx=AnimalIDIdx(~isnan(AnimalIDIdx));
secIDX=length(AnimalIDIdx);
AnimalIDIdx=AnimalIDIdx(1:secIDX);
% figure
TCia=cell(length(AnimalIDIdx), 1);
nearNeighsia=cell(length(AnimalIDIdx), 1);
edgeDistSum=[0,0];
if options.xyNoiseStd>0

    if options.isNoiseCTS

        cellTAll=DTable.Subtype;
        notFixedTypes=~ismember(cellTAll, options.fixedTypes);

        cellTAll(notFixedTypes)=cellTAll(notFixedTypes)+rand(sum(notFixedTypes),1);
        [~, iu, ju]=unique(cellTAll);
        xyNoiseU=options.xyNoiseStd*(rand(length(iu), 2)-0.5);
        xyNoise=xyNoiseU(ju, :);


    else
        xyNoise=options.xyNoiseStd*(rand(size(DTable, 1), 2)-0.5);
    end
else
    xyNoise=zeros(size(DTable, 1), 2);
end

xshiftStep=15000;
yshiftStep=10000;
numNodes=0;
for iAnimalID=1:length(AnimalIDIdx)
    AnimalIDi=AnimalIDIdx(iAnimalID);
    loopFlafID=(AnimalID==AnimalIDi);
    T1=DTable(loopFlafID, :);
    xyNoise1=xyNoise(loopFlafID, :);

    sectionid_codes=T1.SectionNumber;
    sectionid_codesU=unique(sectionid_codes);
    yshift=(AnimalIDi-1)*yshiftStep;
    [~, sectionid_codessInAll]=ismember(sectionid_codesU, sectionU);
    TCib=cell(length(sectionid_codesU), 1);

    % for isectionid_codes=1:min(1,length(sectionid_codesU))
    sectCnt=0;
    nearNeighsib=cell(length(sectionid_codesU), 1);
    for isectionid_codes=1:length(sectionid_codesU)

        

       

        xshift=sectCnt*xshiftStep;
        loopFlagBrg=(sectionid_codes==sectionid_codesU(isectionid_codes));
        T2=T1(loopFlagBrg, :);
        xyNoise2=xyNoise1(loopFlagBrg, :);

        % randNoiseIdx=randi(size(xyNoise, 1), size(T2, 1), 1);

        xcoords=T2.Remapped_X+xyNoise2(:, 1);
        xcoords=xcoords-mean(xcoords);

        ycoords=T2.Remapped_Y+xyNoise2(:, 2);
        ycoords=ycoords-mean(ycoords);
        % xyNoise(1:sum(loopFlagBrg), :)=[];

        % Neuron_cluster_ID


        subtypeVec=T2.Subtype;

        cellSectionVec=ones(length(subtypeVec),1);

        cellSectionVec(:)=(AnimalIDi-1)*length(sectionU)+sectionid_codessInAll(isectionid_codes);

        %
        %         xcoords=xcoords(1:5);
        %         ycoords=ycoords(1:5);
        %         subtype_codesVec=subtype_codesVec(1:5);

        xcoordsTotal=[xcoordsTotal;xcoords+xshift];
        ycoordsTotal=[ycoordsTotal;ycoords+yshift];
        subtypeVecTotal=[subtypeVecTotal;subtypeVec];
        cellSectionVecTotal=[cellSectionVecTotal;cellSectionVec];

        % xcoords=xcoords(subtype_codesVec~=15);
        % ycoords=ycoords(subtype_codesVec~=15);
        % subtype_codesVec=subtype_codesVec(subtype_codesVec~=15);

        % %

        numcells=length(subtypeVec);


        zcoords=[xcoords, ycoords];
        

        subtypeVecU=unique(subtypeVec);
        % figure
        % hold on
        % for ii=1:length(subtype_codesVecU)
        %     zcoordsi=zcoords(subtype_codesVec==subtype_codesVecU(ii), :);
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

        % figure
        % triplot(Delauney_Triangle)

        rEps=options.rEpsilon;
        neighsi= rangesearch(zcoords,zcoords, rEps);
    
        if numNodes>0
            for iNd=1:length(neighsi)
                neighsi{iNd}=neighsi{iNd}+numNodes;
            end
        end
        numNodesi=length(xcoords);
        numNodes=numNodes+numNodesi;
           
        nearNeighsib{isectionid_codes}=neighsi;


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
        TCib{isectionid_codes}=T2;


    sectCnt=sectCnt+1;

    end

    TCia{iAnimalID} = vertcat(TCib{:});
    nearNeighsia{iAnimalID}=vertcat(nearNeighsib{:});
end

TOut=vertcat(TCia{:});

TOut.Remapped_X=xcoordsTotal;
TOut.Remapped_Y=ycoordsTotal;


subtypeVec=subtypeVecTotal;
if isUniformWeight
    Weigths=ones(size(edgesUniqueTotal, 1),1);
else
    Weigths=WeigthsTotal;
end

nearNeighs=vertcat(nearNeighsia{:});

G=graph(edgesUniqueTotal(:,1),edgesUniqueTotal(:,2), Weigths);

G.Nodes.label=[subtypeVec, cellSectionVecTotal];
G.Nodes.Coordinates=([xcoordsTotal, ycoordsTotal]);
if ~isempty(nearNeighs)
    G.Nodes.nearNeighs=nearNeighs;
end

subtypeVecU=unique(subtypeVec);

gStruct=getTransitionStruct(G);
% gStruct.coordinates=G.Nodes.Coordinates;
% gStruct.labels=G.Nodes.label;


nodesAll=(1:length(subtypeVec));

edgeDistStd=edgeDistSum(1)/edgeDistSum(2);

if plotGraph
    figure('WindowState' ,'maximize')
    ncolor=colormap(jet(length(subtypeVecU)));

    h=plot(G,'XData',xcoordsTotal,'YData',ycoordsTotal);
    for iu=1:length(subtypeVecU)
        % ncolor=rand(1, 3);
        if subtypeVecU(iu)==6
            highlight(h,nodesAll(subtypeVec==subtypeVecU(iu)),'NodeColor','k')
        else
            highlight(h,nodesAll(subtypeVec==subtypeVecU(iu)),'NodeColor',ncolor(iu, :))
        end
    end
        hold on
        ax=zeros(length(subtypeVecU), 1);

    for iu=1:length(subtypeVecU)
        ax(iu)=plot(NaN,NaN,'.', 'MarkerFaceColor', ncolor(iu, :), 'MarkerEdgeColor', ncolor(iu, :), 'markersize', 20); %plotting invisible points of desired colors
    end 

    legendText=num2str((1:length(subtypeVecU)).');

    % legend(ax,legendText,'Orientation','horizontal')
    legend(ax,legendText,'NumColumns',4, 'Location', 'Best')


    set(gca,'xtick',sectionU*xshiftStep,'xticklabel',sectionU) 
    set(gca,'ytick',(AnimalIDIdx-1)*yshiftStep,'yticklabel',AnimalIDIdx)
    xlabel('sections')
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
        %     highlight(h,nodesAll(subtype_codesVec==fixedTypes(iFx)),'NodeColor','g')
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
        %     set(gca,'xtick',sectionid_codesAllU*1e5/2,'xticklabel',sectionid_codesAllU) 
        %     set(gca,'ytick',(AnimalIDIdx(1:secIDX)-1)*3000,'yticklabel',(1:secIDX))
        %     xlabel('sectionid_codes')
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
