function [T, imgSize]=readImgQ(imgData,maskData, options)

% resolution=options.resolution;
inmap=options.inmap;

X_no_dither = rgb2ind(imgData,inmap,'nodither');



% imgData=imread(imgData);
imgData=X_no_dither;
% imgData=double(imgData);





imgDataRound=double(X_no_dither)+1;

imgDataRoundV=imgDataRound(:);

mask_dilated = imdilate(maskData>0,ones(5,5));

mask_eroded = imerode(maskData>0,ones(5,5));


BW0 = edge(maskData>0);
BWD = edge(mask_dilated>0);
BWE = edge(mask_eroded>0);

% figure
% subplot(131)
% imshow(BW0)
% subplot(132)
% imshow(BWD)
% subplot(133)
% imshow(BWE)

maskVFlag=BWD(:)|BWE(:);

% maskVFlag=maskVFlag>0;


[~, ~, jU]=unique(imgDataRoundV);


dxy=size(imgData);
nx=dxy(2);
ny=dxy(1);

center_x=repelem((1:nx),ny, 1);
center_x=center_x(:);

center_y=-repmat((1:ny).',nx, 1);

jU=jU(maskVFlag);
center_x=center_x(maskVFlag);
center_y=center_y(maskVFlag);


T=table(center_x, center_y, jU,repelem(options.sampleName, length(jU), 1),'VariableNames', {'Centroid_X','Centroid_Y','cellType','SID'});


imgSize=dxy(1:2);


