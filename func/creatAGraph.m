
function [G,gStruct,  h, BIDAllU, AIDIdx, TOut, edgeDistStd]=creatAGraph(filename, options)

T = readtable(filename);

edgesUniqueTotal=[];
numcellsTotal=0;
xcoordsTotal=[];
ycoordsTotal=[];
cellSubtypeVecTotal=[];
cellSectionVecTotal=[];
AID=T.AID;
BIDAll=T.BID;
BIDAllU=unique(BIDAll);
AIDIdx=(unique(AID)).';
AIDIdx=AIDIdx(~isnan(AIDIdx));
TCia=cell(length(AIDIdx), 1);
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

nearNeighsia=cell(length(AIDIdx),1);

for iAID=1:length(AIDIdx)
    AIDi=AIDIdx(iAID);
    loopFlafID=(AID==AIDi);
    T1=T(loopFlafID, :);
    xyNoise1=xyNoise(loopFlafID, :);

    BID=T1.BID;
    BIDU=unique(BID);
    if length(BIDU)<2
        yshift=(iAID-1)*3000;
    else
        yshift=(AIDi-1)*3000;
    end
    [~, BInAll]=ismember(BIDU, BIDAllU);
    TCib=cell(length(BIDU), 1);
    nearNeighsib=cell(length(BIDU), 1);

    for iBregma=1:length(BIDU)       

        xshift=(BIDU(iBregma))*1e5/2;
        loopFlagBrg=(BID==BIDU(iBregma));
        T2=T1(loopFlagBrg, :);
        xyNoise2=xyNoise1(loopFlagBrg, :);
        maxX=max(T2.Centroid_X);
        minX=min(T2.Centroid_X);     
        xcoords=T2.Centroid_X+xyNoise2(:, 1);
        gFlag=xcoords>maxX;
        lFlag=xcoords<minX;    
        xcoords(gFlag|lFlag)=xcoords(gFlag|lFlag)-0.9*xyNoise2(gFlag|lFlag, 1);
        xcoords=xcoords-mean(xcoords);
        maxY=max(T2.Centroid_Y);
        minY=min(T2.Centroid_Y);
        ycoords=T2.Centroid_Y+xyNoise2(:, 2);
        gFlag=ycoords>maxY;
        lFlag=ycoords<minY;
        ycoords(gFlag|lFlag)=ycoords(gFlag|lFlag)-0.9*xyNoise2(gFlag|lFlag, 2);
        ycoords=ycoords-mean(ycoords);
        cellSubtypeVec=T2.cellSubtype;

        cellSectionVec=ones(length(cellSubtypeVec),1);

        cellSectionVec(:)=(AIDi-1)*length(BIDAllU)+BInAll(iBregma);
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
            if options.rEps>0
                rEps=options.rEps;
            else
                numNeighs=2;
                [~, D] = knnsearch(zcoords,zcoords, 'K',numNeighs);
                DS=sort(D(:));
                rEps=DS(end-10);
            end
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

        edges=sort(edges, 2);
        [edgesUnique,~,~]=unique(edges, 'rows', 'stable');

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
        TCib{iBregma}=T2;

    end

    TCia{iAID} = vertcat(TCib{:});
    nearNeighsia{iAID}=vertcat(nearNeighsib{:});
end

TOut=vertcat(TCia{:});

TOut.Centroid_X=xcoordsTotal;
TOut.Centroid_Y=ycoordsTotal;


cellSubtypeVec=cellSubtypeVecTotal;
Weigths=ones(size(edgesUniqueTotal, 1),1);


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

edgeDistStd=edgeDistSum(1)/edgeDistSum(2);

