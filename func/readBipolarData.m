
function To=readBipolarData(folderName, options)

numFileSection=options.retinaAndSection;


if options.isRandom


    numCTypesRand=options.nCTypes;
    numNodesRand=options.nNodesRand;
    numNodesHold=floor(numNodesRand/5);

    sectNumberIdx=[1, 2];
    Tc=cell(2, 1);
    for sectNumberi=sectNumberIdx
        yshift=(sectNumberi-1)*sqrt(numNodesRand);

        if sectNumberi==1
            numNodesRandi=numNodesHold;
        else
            numNodesRandi=numNodesRand;
        end

        xcoords=sqrt(numNodesRandi)*rand(numNodesRandi,1);
        xcoords=xcoords-mean(xcoords);
        ycoords=sqrt(numNodesRandi)*rand(numNodesRandi,1);
        ycoords=ycoords-mean(ycoords);

        cellSubtypeVec=randsample(numCTypesRand,numNodesRandi, true);
        cellSectionVec=ones(length(cellSubtypeVec),1);

        cellSectionVec(:)=sectNumberi;

        Ti.Remapped_X=xcoords;

        Ti.Remapped_Y=ycoords+yshift;

        Ti.Subtype=cellSubtypeVec;
        Ti.SectionNumber=cellSectionVec;
        Tc{sectNumberi}=Ti;


    end
    To=vertcat(Tc{:});

else


    TfC=cell(numFileSection(1), 1);

    for ifile=1:numFileSection(1)
        csvName=strcat('Retina', num2str(ifile));
        filename=strcat(folderName,'\',csvName, '.csv');
        % filename='..\data\Retina1.csv';

        Tfi = readtable(filename);
        sectNumber=Tfi.SectionNumber;

        % cellSubtypeVec=cellSubtypeVec(randperm(length(cellSubtypeVec)));

        % iterMax=100;
        % s_cellSubtypeVec=zeros(numcells, iterMax);
        % for ii =1:iterMax
        %     s_cellSubtypeVec(:,ii)=cellSubtypeVec(randperm(numcells));
        % end

        Tfi=Tfi(~isnan(sectNumber), :);

        Tfi.animalID=repelem(ifile,size(Tfi,1), 1);
        TfC{ifile}=Tfi;

    end

    To=vertcat(TfC{:});

end
