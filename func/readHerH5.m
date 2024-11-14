function [allData, structNames]=readHerH5(filename)

info = h5info(filename);

namesOut=readDSNames(info, []);

GIdxs=find(namesOut=='$');
numGs=length(GIdxs);
allGroups=cell(numGs, 1);
GIdxs0=1;
for iGId=1:length(GIdxs)

    allGroupsi=namesOut(GIdxs0:GIdxs(iGId)-1);
    DIdxs=find(allGroupsi=='&');

    if ~isempty(DIdxs)
        DIdxsE=[DIdxs, length(allGroupsi)+1];
        if length(DIdxs)==1
            gnamei=allGroupsi(1:DIdxsE(2)-1);
            gnamei(gnamei=='&')=[];
            DIdxs0=DIdxsE(2)+1;
        else
            gnamei=allGroupsi(1:DIdxsE(1)-1);
            DIdxs0=DIdxsE(1)+1;
        end           

        numDSi=length(DIdxs);
        allDatasetsi=cell(1, numDSi);
        
        for iDId=1:numDSi
            allDatasetsi{iDId}=strcat(gnamei, allGroupsi(DIdxs0:DIdxsE(iDId+1)-1));
            DIdxs0=DIdxsE(iDId+1)+1;
        end
            allGroups{iGId}=allDatasetsi;
    else
        allGroups{iGId}='';
    end






    GIdxs0=GIdxs(iGId)+1;
end

allNames=cat(2, allGroups{:});
numData=length(allNames);

allData=cell(numData, 1);


dataTable = table;

structNames=cell(numData, 1);

for ii =1:numData
    namei=allNames{ii};
    namei(namei=='/')='_';
    namei=namei(2:end);
    structNames{ii}=namei;
    allDatai.name=namei;
    datai=h5read(filename,allNames{ii});
    allDatai.data=datai;
    allData{ii}=allDatai;

    nameif=find(namei=='_', 1);
    if ~isempty(nameif)
        namei=namei(nameif+1:end);
    end

    if size(datai, 1)==114188
        dataTable.(namei)=double(datai);
    elseif size(datai, 2)==114188
         dataTable.(namei)=double(datai.');
    % else
    %     dataTable.(namei)=double(datai);
    end
end

check=1;


