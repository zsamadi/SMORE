function T=readCBRTable(fileName, options)

% cd(fileparts(which(mfilename)));
% cd('..\')
%
% addpath(genpath(pwd))
%
% dataFolder='D:\nucla\P1\matIO\data\BRACS\BRACS_RoI\latest_version\test\0_N\';
%
% fileName=strcat(dataFolder, 'BRACS_1286_N_40_labels_props.csv');
%

T = readtable(fileName);
clusterIdx=T.cluster;

if options.isHighlight
    outputFigName=options.outputFigName;
    imgName=strcat(fileName(1:end-10), '.tif');

    labels=T.label;

    thisImage = imread(imgName);


    % clusterColors=[0 0.4470 0.7410;0.9290 0.6940 0.1250;0.8500 0.3250 0.0980;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;0 0 0;0 1 1;1 0 1;0 1 0;1 0.64 0;1 0 0;rand(numClusters, 3)];

    clusterColors=options.clusterColors;


    coloredImage=ones([size(thisImage), 3]);



    for il=1:length(labels)

        [row, column]=find(thisImage==il);
        clusterColorsi=clusterColors(clusterIdx(il), :);

        for irc=[row, column].'


            % aa=repmat(clusterColors(clusterIdx(il), :), length(row)* length(row), 1);
            % aa1=reshape(aa, length(row), length(row), 3);

            coloredImage(irc(1), irc(2), :)=clusterColorsi;
        end
    end

    figure
    imshow(coloredImage)

    figname=outputFigName;

    
    saveas(gcf,figname);



    thisImageBW=double(thisImage>0);


    figure
    imshow(thisImageBW)
end



check=1;















