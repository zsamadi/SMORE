function ctShuffle=getShuffleTypes(cTypes, config)

cTypeSize=size(cTypes);
cTypes=cTypes(:);
fixedTypes=config.fixedTypes;
if ~isempty(fixedTypes)
    fixedIndex=ismember(cTypes, fixedTypes);
else
    fixedIndex=false(size(cTypes));
end

fixedIndex(config.fixedNodes)=true;
ctShuffleC=cell(config.numShuffle, 1);

if strcmpi(config.shuffleMode, "kernel")

    nearNeighs=config.nearNeighs;

    for iShuff=1:config.numShuffle
        ctShuffle=cTypes;
        nodesAll=(1:length(cTypes));
        nodesAllNotFixed=nodesAll(~fixedIndex);

        for iGN=nodesAllNotFixed
            neighsIGN=nearNeighs{iGN};
            neighsIGNE=[iGN;neighsIGN(:)];
            cPTypesNgi=cTypes(neighsIGNE);

            cPTypesNgi(ismember(cPTypesNgi, fixedTypes))=[];


            randIndex=randi(length(cPTypesNgi));

            ctShuffle(iGN)=cPTypesNgi(randIndex);
        end

        ctShuffleC{iShuff}=ctShuffle;

    end

    ctShuffle=vertcat(ctShuffleC{:});

elseif strcmpi(config.shuffleMode, "kernelPath")


    nearNeighs=config.nearNeighs;

    for iShuff=1:config.numShuffle
        ctShuffle=cTypes;
        nodesAll=(1:length(cTypes));
        nodesAllNotFixed=nodesAll(~fixedIndex);

        for iGN=nodesAllNotFixed
            neighsIGN=nearNeighs{iGN};
            neighsIGNE=neighsIGN(:);
            cPTypesNgi=cTypes(neighsIGNE);










            cPTypesNgi(ismember(cPTypesNgi, fixedTypes))=[];

            if isempty(cPTypesNgi)
                check=1;
            end

             cPTypesi=cPTypesNgi(1);

            if ~isempty(config.partitionInd)
                partitionIndi=config.partitionInd(cPTypesi);
                partitionSeti=config.partitionSet{partitionIndi};
                cPTypesNgi(~ismember(cPTypesNgi, partitionSeti))=[];
            end



            randIndex=randi(length(cPTypesNgi));

            ctShuffle(iGN)=cPTypesNgi(randIndex);
        end

        ctShuffleC{iShuff}=ctShuffle;

    end

    ctShuffle=vertcat(ctShuffleC{:});


else


    isSectShuffle=config.isSectShuffle;
    cSections=config.cSections;
    cSections=cSections(:);
    % xyCoordinates=config.xyCoordinates;




    cSections=cSections(~fixedIndex);
    % xyCoordinates=xyCoordinates(~fixedIndex, :);
    ctShuffle0=cTypes(~fixedIndex);

    % kCluster=config.numClusters;

    % cSectionsU=unique(cSections);
    numNodes=length(cSections);
    cSections2=ones(numNodes, 2);


    cSections2(:, 1)=cSections;



    cSections2U=unique(cSections2, 'rows');

    for iShuff=1:config.numShuffle

        ctShuffle=cTypes;


        if isSectShuffle

            for iSec=1:size(cSections2U, 1)
                secFlag=all(cSections2==cSections2U(iSec, :),2);

                ctShufflet=ctShuffle0(secFlag);

                ctShufflet=ctShufflet(randperm(length(ctShufflet)));
                ctShuffle0(secFlag)=ctShufflet;
            end
        else

            ctShuffle0=ctShuffle0(randperm(length(ctShuffle0)));
        end

        ctShuffle(~fixedIndex)=ctShuffle0;
        ctShuffleC{iShuff}=ctShuffle;
    end
    ctShuffle=vertcat(ctShuffleC{:});

end

if cTypeSize(1)<cTypeSize(2)
    ctShuffle=ctShuffle.';
end




