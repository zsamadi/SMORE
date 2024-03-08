
function [G,gStruct, TOut, nodeTypeVecU, haveZ]=creatAGraph(filename, options)

T = readtable(filename);

if options.ND==0
    varNames=T.Properties.VariableNames;
    if strcmpi(varNames{5}, 'Centroid_Z')
        haveZ=true;
    else
         haveZ=false;
    end
elseif  options.ND==3
    haveZ=true;
else
 haveZ=false;
end




AID=T.NID;
BIDAll=ones(size(AID));
BIDAllU=unique(BIDAll);
AIDIdx=(unique(AID)).';
AIDIdx=AIDIdx(~isnan(AIDIdx));
TCia=cell(length(AIDIdx), 1);
edgeDistSum=[0,0];
numNodes=0;

nearNeighsia=cell(length(AIDIdx),1);


edgesUniqueTotal=cell(length(AIDIdx),1);
numcellsTotal=0;
coordsTotal=cell(1000,1);
nodeTypeVecTotal=cell(1000,1);
cellSectionVecTotal=cell(1000,1);


iABID=0;
xshift=0;
yshift=0;

for iAID=1:length(AIDIdx)
    AIDi=AIDIdx(iAID);
    loopFlafID=(AID==AIDi);
    T1=T(loopFlafID, :);

    % BID=T1.BID;
    BID=ones(sum(loopFlafID), 1);

    BIDU=unique(BID);

    [~, BInAll]=ismember(BIDU, BIDAllU);
    TCib=cell(length(BIDU), 1);
    nearNeighsib=cell(length(BIDU), 1);
    xshift=0;

    for iBID=1:length(BIDU)  
        iABID=iABID+1;

        loopFlagBID=(BID==BIDU(iBID));
        T2=T1(loopFlagBID, :);
        xcoords=T2.Centroid_X; 
        xcoords=xcoords-mean(xcoords);
        ycoords=T2.Centroid_Y;
        ycoords=ycoords-mean(ycoords);

        if haveZ
            zcoords=T2.Centroid_Z;
            zcoords=zcoords-mean(zcoords);
        end




        nodeTypeVec=T2.nodeType;

        cellSectionVec=ones(length(nodeTypeVec),1);

        cellSectionVec(:)=(AIDi-1)*length(BIDAllU)+BInAll(iBID);
        
        nodeTypeVecTotal{iABID}=nodeTypeVec;
        cellSectionVecTotal{iABID}=cellSectionVec;
        numcells=length(nodeTypeVec);
        if haveZ
            xyzcoords=[xcoords, ycoords, zcoords];
            coordsTotal{iABID}=[xcoords+xshift, ycoords+yshift, zcoords];

        else
            xyzcoords=[xcoords, ycoords];
            coordsTotal{iABID}=[xcoords+xshift, ycoords+yshift];

        end

        if strcmpi(options.gMode,"delaunay")


            Delauney_Triangle=delaunayTriangulation(xyzcoords);

            connectTriangle=Delauney_Triangle.ConnectivityList;

            edges=[connectTriangle(:,1:2); connectTriangle(:,2:3); connectTriangle(:,[1,3])];
        elseif strcmpi(options.gMode,"knn")

            numNeighs=options.numNeighs;

            nearNeighs = knnsearch(xyzcoords,xyzcoords, 'K',numNeighs);

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
                [~, D] = knnsearch(xyzcoords,xyzcoords, 'K',numNeighs);
               
                rEps=4*mean(D(D>0));
            end
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

           rEps=options.rEpsilon;
           neighsi= rangesearch(xyzcoords,xyzcoords, rEps);
    
           if numNodes>0
                for iNd=1:length(neighsi)
                     neighsi{iNd}=neighsi{iNd}+numNodes;
                end
           end
           numNodesi=length(xcoords);
           numNodes=numNodes+numNodesi;
           
           nearNeighsib{iBID}=neighsi;

        edges=sort(edges, 2);
        [edgesUnique,~,~]=unique(edges, 'rows', 'stable');

        edgeDist02=xyzcoords(edgesUnique(:,2), :)-xyzcoords(edgesUnique(:,1), :);
        edgeDist02=edgeDist02.^2;
        edgeDist02=sum(edgeDist02,2);
        Weigths=exp(-edgeDist02/mean(sqrt(edgeDist02))^2);     

        weightThreshold=1e-10;
        edgeDistSum=edgeDistSum+[sum(sqrt(edgeDist02(Weigths>weightThreshold))), sum(Weigths>weightThreshold)];
        edgesUnique=edgesUnique(Weigths>weightThreshold, :);
        edgesUniqueTotal{iABID}=edgesUnique+numcellsTotal;
        numcellsTotal=numcellsTotal+numcells;
        TCib{iBID}=T2;
        xshift=xshift+2.2*max(abs(coordsTotal{iABID}(:, 1)));


    end
    yshift=yshift+2.2*max(abs(coordsTotal{iABID}(:, 2)));

    TCia{iAID} = vertcat(TCib{:});
    nearNeighsia{iAID}=vertcat(nearNeighsib{:});
end


coordsTotal=vertcat(coordsTotal{:});
nodeTypeVecTotal=vertcat(nodeTypeVecTotal{:});
cellSectionVecTotal=vertcat(cellSectionVecTotal{:});

edgesUniqueTotal=vertcat(edgesUniqueTotal{:});

TOut=vertcat(TCia{:});

TOut.Centroid_X=coordsTotal(:, 1);
TOut.Centroid_Y=coordsTotal(:, 2);

if haveZ
    TOut.Centroid_Z=coordsTotal(:, 3);
end
    

nodeTypeVec=nodeTypeVecTotal;
Weigths=ones(size(edgesUniqueTotal, 1),1);


% nearNeighs=vertcat(nearNeighsia{:});
nearNeighs=[];

G=graph(edgesUniqueTotal(:,1),edgesUniqueTotal(:,2), Weigths);

G.Nodes.label=[nodeTypeVec, cellSectionVecTotal];
G.Nodes.Coordinates=coordsTotal;
if ~isempty(nearNeighs)
    G.Nodes.nearNeighs=nearNeighs;
end

[nodeTypeVecU, ~, jU]=unique(nodeTypeVec);
cellCnts=accumarray(jU, 1);
[~, ist]=sort(cellCnts, 'descend');
nodeTypeVecU=nodeTypeVecU(ist);
gStruct=getTransitionStruct(G);

