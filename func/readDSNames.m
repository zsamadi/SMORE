function namesOut=readDSNames(DSet, namesIn)
namesOut=namesIn;
dset=DSet.Datasets;
if ~isempty(dset)
    for id=1:length(dset)
        dseti=dset(id);
        namesOuti=dseti.Name;
        namesOut=strcat(namesOut, '&/', namesOuti);
    end
end
groupID=DSet.Groups;
if ~isempty(groupID)
    for ig=1:length(groupID)
        namesOuti=readDSNames(groupID(ig), groupID(ig).Name);
        namesOut=strcat(namesOut, '$', namesOuti);
    end
end




