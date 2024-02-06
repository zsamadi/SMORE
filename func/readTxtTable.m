function [T, cellTypesOut, txtSize]=readTxtTable(fileName, noiseVarParam)


txtData=importdata(fileName);
txtData=horzcat(txtData{:});
txDataAll='';
cellTypesOut=char(65:90);
cellTypesOut=cellTypesOut(:);

lineLen=zeros(length(txtData), 1);

    txtDatai=upper(txtData);

    txtDataiFlag=find(txtDatai== ' ');

lineMaxLen=floor(sqrt(length(txtData)));



    islined=length(txtDataiFlag)<lineMaxLen;
    txtDataij=cell(ceil(length(txtData)/lineMaxLen),1);
    j=1;

    while ~islined

        nextIDX=find(txtDataiFlag>lineMaxLen, 1);
        txtDataij{j}=txtDatai(1:txtDataiFlag(nextIDX));
        txtDatai(1:txtDataiFlag(nextIDX))=[];
        txtDataiFlag=find(txtDatai== ' ');
        islined=all(txtDataiFlag<=lineMaxLen);
        j=j+1;
    end
    txtDataij{j}=txtDatai;
    txtData=txtDataij(1:j);



for iTD=1:length(txtData)
    txtDatai=txtData{iTD};

    txtDataiFlag=ismember(txtDatai, cellTypesOut);
    txtDatai(~txtDataiFlag)=[];

    txDataAll=strcat(txDataAll,txtDatai);
    lineLen(iTD)=length(txtDatai);
end

txDataAll=txDataAll(:);


center_y=repelem((1:length(lineLen)).', lineLen);
center_y_noise=noiseVarParam*rand(size(center_y));

center_y_noise=center_y_noise-mean(center_y_noise);

center_y=center_y+center_y_noise;

lineBreaks=cumsum(lineLen);
lineBreaks=[0;lineBreaks(1:end-1)];
center_x=(1:length(center_y)).';

center_x_noise=noiseVarParam*rand(size(center_x));

center_x_noise=center_x_noise-mean(center_x_noise);

lineBreaksRep=repelem(lineBreaks, lineLen);

center_x=center_x - lineBreaksRep;

center_x=center_x+center_x_noise;

T  = cell2table(cell(length(center_y),1), 'VariableNames', {'center_x'});


T.Centroid_X=center_x;
T.Centroid_Y=center_y;


T.cellSubtype=double(txDataAll)-64;
txtSize=length(center_x);


