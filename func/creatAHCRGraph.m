
function [G,gStruct, cidOut, cellTypeVecU, haveZ, SIDIdxS, TOut]=creatAHCRGraph(filename, options)

T = readtable(filename);

varNames=T.Properties.VariableNames;

haveZ=any(matches(varNames, 'Centroid_Z'));
[SIDIdxS, ~, AID]=unique(T.SID);
T.AIDO=AID;

if haveZ && ~options.ND

    T.SID=[AID, T.Centroid_Z];
    [~, ~, AID]=unique(T.SID, 'rows');
    haveZ=false;
end




T.CIDNm=(1:size(T, 1)).';



AIDIdx=(unique(AID)).';
AIDIdx=AIDIdx(~isnan(AIDIdx));

TCia=cell(length(AIDIdx), 1);
nearNeighsia=cell(length(AIDIdx),1);
edgesUniqueTotal=cell(length(AIDIdx),1);


NDxy=ceil(sqrt(length(AIDIdx)));




% BIDAll=mod(AID-1, NDxy)+1;
BIDAll=mod(AID-1, NDxy)+1;


AID=ceil(AID/NDxy);

AIDIdx=(unique(AID)).';
AIDIdx=AIDIdx(~isnan(AIDIdx));

BIDAllU=unique(BIDAll);




edgeDistSum=[0,0];
numNodes=0;


numcellsTotal=0;
coordsTotal=cell(1000,1);
cellTypeVecTotal=cell(1000,1);
cellSectionVecTotal=cell(1000,1);


iABID=0;
yshift=0;
numNeighs=options.nNeighs;

for iAID=1:length(AIDIdx)
    AIDi=AIDIdx(iAID);
    loopFlafID=(AID==AIDi);
    T1=T(loopFlafID, :);

    BID=BIDAll(loopFlafID, :);
    % BID=ones(sum(loopFlafID), 1);

    BIDU=unique(BID);

    [~, BInAll]=ismember(BIDU, BIDAllU);
    TCib=cell(length(BIDU), 1);
    nearNeighsib=cell(length(BIDU), 1);
    xshift=0;

    for iBID=1:length(BIDU)  
        iABID=iABID+1;

        loopFlagBID=(BID==BIDU(iBID));
        T2=T1(loopFlagBID, :);
        if haveZ
            [~, iu]=unique([T2.Centroid_X, T2.Centroid_Y, T2.Centroid_Z], 'rows', 'stable');
            T2=T2(iu, :);
        else
            [~, iu]=unique([T2.Centroid_X, T2.Centroid_Y], 'rows', 'stable');
            T2=T2(iu, :);
        end


        xcoords=T2.Centroid_X; 
        xcoords=xcoords-min(xcoords);
        ycoords=T2.Centroid_Y;
        ycoords=ycoords-min(ycoords);

        if haveZ
            zcoords=T2.Centroid_Z;
            zcoords=zcoords-mean(zcoords);

        end




        cellTypeVec=T2.cellType;

        % cellSectionVec=ones(length(cellTypeVec),1);

        % cellSectionVec(:)=(AIDi-1)*length(BIDAllU)+BInAll(iBID);
        cellSectionVec=T2.AIDO;
        
        cellTypeVecTotal{iABID}=cellTypeVec;
        cellSectionVecTotal{iABID}=cellSectionVec;
        numcells=length(cellTypeVec);
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
        TCib{iBID}=T2.CIDNm;
        xshift=1.1*max(abs(coordsTotal{iABID}(:, 1)));
        % figure
        % 
        % scatter(coordsTotal{iABID}(:, 1),coordsTotal{iABID}(:, 2))
        % check=1;


    end
    coordSoFar=vertcat(coordsTotal{:});
    yshift=max(abs(coordSoFar(:, 2)));

    TCia{iAID} = vertcat(TCib{:});
    nearNeighsia{iAID}=vertcat(nearNeighsib{:});
end


coordsTotal=vertcat(coordsTotal{:});
cellTypeVecTotal=vertcat(cellTypeVecTotal{:});
cellSectionVecTotal=vertcat(cellSectionVecTotal{:});

edgesUniqueTotal=vertcat(edgesUniqueTotal{:});

cidOut=vertcat(TCia{:});

% cidOut=cidOut.CID;

% cidOut.Centroid_X=coordsTotal(:, 1);
% cidOut.Centroid_Y=coordsTotal(:, 2);
% 
% if haveZ
%     cidOut.Centroid_Z=coordsTotal(:, 3);
% end
    

cellTypeVec=cellTypeVecTotal;
Weigths=ones(size(edgesUniqueTotal, 1),1);


% nearNeighs=vertcat(nearNeighsia{:});
nearNeighs=[];

G=graph(edgesUniqueTotal(:,1),edgesUniqueTotal(:,2), Weigths);

if size(G.Nodes, 1)<length(cellTypeVec)
    G=addnode(G, length(cellTypeVec)-size(G.Nodes, 1));
end

[cellTypeVecU, ~, jcpU]=unique(cellTypeVec);
% cellCnts=accumarray(jcpU, 1);
% [~, ist]=sort(cellCnts, 'descend');
% cellTypeVecU=cellTypeVecU(ist);

G.Nodes.label=[jcpU, cellSectionVecTotal];
% cidOut.cellType=jcpU;
G.Nodes.Coordinates=coordsTotal;
if ~isempty(nearNeighs)
    G.Nodes.nearNeighs=nearNeighs;
end
TOut=T(cidOut, :);

gStruct=getTransitionStruct(G);

