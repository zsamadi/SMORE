function result=writeCBRTable(folderName, options)
numClusters0=options.numClusters0;
numClusters1=options.numClusters1;


warning off
startPath=folderName;
topLevelFolder = startPath;
% filePattern = sprintf('%s/output/hypoW4_*/*.mat', topLevelFolder);
filePattern = sprintf('%s/**/*.csv', topLevelFolder);

allFileInfo = dir(filePattern);

isFolder = [allFileInfo.isdir]; % Logical list of what item is a folder or not.
% Now set those folder entries to null, essentially deleting/removing them from the list.
allFileInfo(isFolder) = [];
% Get a cell array of strings.  We don't really use it.  I'm just showing you how to get it in case you want it.
listOfFolderNames = unique({allFileInfo.folder});
numberOfFolders = length(listOfFolderNames);
fprintf('The total number of folders to look in is %d.\n', numberOfFolders);

% Get a cell array of base filename strings.  We don't really use it.  I'm just showing you how to get it in case you want it.
listOfFileNames = {allFileInfo.name};
totalNumberOfFiles = length(listOfFileNames);
fprintf('The total number of files in those %d folders is %d.\n', numberOfFolders, totalNumberOfFiles);

% Process all files in those folders.
totalNumberOfFiles = length(allFileInfo);
% Now we have a list of all files, matching the pattern, in the top level folder and its subfolders.
TAll=cell(totalNumberOfFiles, 1);
TLen=zeros(totalNumberOfFiles, 1);
areaMeanV=zeros(totalNumberOfFiles, 1);
parfor k = 1 : totalNumberOfFiles

    % Go through all those files.
    thisFolder = allFileInfo(k).folder;
    thisBaseFileName = allFileInfo(k).name;
    fullFileName = fullfile(thisFolder, thisBaseFileName);

    T = readtable(fullFileName);

    dashPos=find(thisFolder=='_');

    TLenk=size(T, 1);


        T.cancerType=string(repelem(thisFolder(dashPos(end)+1:end), TLenk, 1));
        T.sampleID=repelem(k, TLenk, 1);


    TAll{k}=T;
    TLen(k)=TLenk;

    areaMeanV(k)=mean(T.area)
end

TAll=vertcat(TAll{:});

%% cluster

TAllFileName=strcat('D:\nucla\P1\matIO\data\', 'cellTable.csv');







% featureMatrix=[TAll.axis_major_length, TAll.axis_minor_length,TAll.area, TAll.solidity,   TAll.eccentricity,TAll.perimeter,TAll.equivalent_diameter_area, TAll.extent, TAll.feret_diameter_max,TAll.orientation];

intensity_min=mean([TAll.intensity_min_0, TAll.intensity_min_1, TAll.intensity_min_2], 2);

intensity_max=mean([TAll.intensity_max_0, TAll.intensity_max_1, TAll.intensity_max_2], 2);

intensity_mean=mean([TAll.intensity_mean_0, TAll.intensity_mean_1, TAll.intensity_mean_2], 2);


eccentAll=TAll.eccentricity;

QVals = quantile(eccentAll,[1/3, 2/3]);


eccentAllQ=double((eccentAll>QVals(2)))+2;


eccentAllQ(eccentAll<QVals(1))=1;




% eccentAll=eccentAll.^1.5;

eccentAll=eccentAll-min(eccentAll);


eccentAll=eccentAll/(max(eccentAll)+eps);

eccentAll=eccentAll.^2;



% eccentAllQ=round((eccentAll+0.5/numClusters0)*numClusters0);



areaAll=TAll.area;

areaMeanVE=repelem(areaMeanV, TLen,1);

areaAll=areaAll./areaMeanVE;

% 
% QVals = quantile(areaAll,[1/3, 2/3]);
% 
% 
% areaAllQ=double((areaAll>QVals(2)))+2;
% 
% 
% areaAllQ(areaAll<QVals(1))=1;






QMax = quantile(areaAll,0.99);
areaAll(areaAll>QMax)=QMax;

QMin = quantile(areaAll,0.01);
areaAll(areaAll<QMin)=QMin;

areaAll=areaAll-QMin;


areaAll=areaAll/(QMax-QMin+eps);



areaAll=sqrt(areaAll);

areaAllQ=round((areaAll+0.5/numClusters1)*numClusters1);




featureMatrix0=[TAll.eccentricity,TAll.area];


% numFeatures=size(featureMatrix0, 2);
% % numFeatures=3;
% 
% featureMatrix0=featureMatrix0(:, 1:numFeatures);
featureMatrix0=featureMatrix0-mean(featureMatrix0);

featureMatrix0=featureMatrix0./std(featureMatrix0);

% clusterIdx=1+double(featureMatrix>0).*2.^(numFeatures-1:-1:0);




clusterIdx0=kmeans(featureMatrix0,numClusters0);



featureMatrix1=[intensity_min,intensity_mean,intensity_max];

featureMatrix1=featureMatrix1-mean(featureMatrix1);

featureMatrix1=featureMatrix1./std(featureMatrix1);





clusterIdx1=kmeans(featureMatrix1,numClusters1);

% clusterIdx01=[clusterIdx0,clusterIdx1];

clusterIdx01=[areaAllQ, eccentAllQ];

[~, ~, clusterIdx01j]=unique(clusterIdx01, 'rows');



clusterIdx01j(clusterIdx01j==7)=8;

clusterIdx01j(clusterIdx01j==2)=1;
clusterIdx01j(clusterIdx01j==3)=1;


TAll.cluster=clusterIdx01j;




%% write output
writetable(TAll(:, end-2:end),TAllFileName);

startIndex=1;

for k = 1 : totalNumberOfFiles

    % Go through all those files.
    thisFolder = allFileInfo(k).folder;
    thisBaseFileName = allFileInfo(k).name;
    fullFileName = fullfile(thisFolder, thisBaseFileName);

    Tk=TAll(startIndex:startIndex+TLen(k)-1, :);

    writetable(Tk,fullFileName);
    startIndex=startIndex+TLen(k);

end




check=1;


result=0;












