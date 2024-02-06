function cTypesKernel=getKernelShuffleTypes(cTypes,config, options)

fixedTypes=options.fixedTypes;

nearNeighs=config.nearNeighs;

cTypesKernel=cTypes;


fixedIndex=ismember(cTypes, fixedTypes);
nodesAll=(1:length(cTypes));
nodesAllNotFixed=nodesAll(~fixedIndex);

for iGN=nodesAllNotFixed
    neighsIGN=nearNeighs{iGN};
    neighsIGNE=[iGN;neighsIGN(:)];
    cPTypesi=cTypes(neighsIGNE);

    cPTypesi(ismember(cPTypesi, fixedTypes))=[];


    randIndex=randi(length(cPTypesi));

    cTypesKernel(iGN)=cPTypesi(randIndex);
end
