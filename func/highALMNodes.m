function highALMNodes(G, ppHMNodesCell,outputFolderName, config)


offsetHL=0;
ndesPlotIDx=config.ndesPlotIDx;
fignamExtnd=config.fignamExtnd;

numHighlights=config.numHLights;
xcoordsTotal=G.Nodes.Coordinates(:, 1);
ycoordsTotal=G.Nodes.Coordinates(:, 2);

if config.is3D
    zcoordsTotal=G.Nodes.Coordinates(:, 3);
end

sidS=G.Nodes.label(:, 2);


[sidSU, ~, jsidS]=unique(sidS);

% xShift=11;
% yShift=6.5;
xShift=2000;
yShift=2000;

numCols=ceil(sqrt(length(sidSU)));

for iU=1:length(sidSU)
    selectI=jsidS==iU;
    xcoordsTotal(selectI)=xcoordsTotal(selectI)-min(xcoordsTotal(selectI))+xShift*(mod(iU-1, numCols));
    ycoordsTotal(selectI)=ycoordsTotal(selectI)-min(ycoordsTotal(selectI))+yShift*floor((iU-1)/numCols);
        % ycoordsTotal(selectI)=ycoordsTotal(selectI)-min(ycoordsTotal(selectI))+yShift*(mod(iU-1,4));

end

    


% numHighlights=length(ppHMNodesCell);

f=figure('visible','on');
f.Position(3:4) = f.Position(3:4)*1.5;
f.Position(1:2) = f.Position(1:2)/1.5;
grColor=[0.8, 0.8, 0.8];

if config.is3D
    hAll=plot(G,'XData',xcoordsTotal,'YData',ycoordsTotal,'ZData',zcoordsTotal,'NodeColor',grColor, 'EdgeColor','none', 'MarkerSize',1);
    % hAll=plot(G,'XData',xcoordsTotal,'YData',ycoordsTotal,'ZData',zcoordsTotal);
else
    hAll=plot(G,'XData',xcoordsTotal,'YData',ycoordsTotal,'NodeColor',grColor,  'EdgeColor','none', 'MarkerSize',1);
end

ncolor = [
    0.0, 1.0, 0.0;    % Green
    1.0, 0.0, 1.0;    % Magenta
    1.0, 0.647, 0.0;  % Orange
    1.0, 0.0, 0.0;    % Red
    0, 0, 0;
    0.0, 1.0, 1.0;    % Cyan
    % 0.0, 0.0, 1.0;    % Blue
   
    0.502, 0.0, 0.502; % Purple
    0.8, 1,0;    % Yellow
    0.0, 0.502, 0.0;  % Dark Green
    0.647, 0.165, 0.165; % Brown
];

% ncolor=[0 1 0;1 0 1;1 0.64 0;1 0 0;0 0 0;0 1 1];
ncolorAll=[0.9290 0.6940 0.1250;0.7500 0.3250 0.0980;0.4940 0.1840 0.5560;0.3010 0.7450 0.9330;0.4660 0.6740 0.1880;0.6350 0.0780 0.1840;rand(60,3)];

ncolorAll=[ncolor;ncolorAll];

% ncolorAll=[1 0 1;0 1 0;1 0.64 0;1 0 0;0 0 0;0 1 1;rand(numHighlights,3)];

for iMG=1:numHighlights

    clMi=ncolorAll(iMG, :);

    motifPaths=ppHMNodesCell{iMG+offsetHL};

    motifPaths=motifPaths(:);
    if ~isempty(ndesPlotIDx)
        [~, motifPaths]=ismember(motifPaths, ndesPlotIDx);
        motifPaths=motifPaths(motifPaths>0);
    end



    % highlight(hAll,motifPaths(:),'NodeColor',clMi, 'EdgeColor',clMi);
    highlight(hAll,motifPaths(:),'NodeColor',clMi);

end

    grid minor

    xlabel('X')
    ylabel('Y')
    numLegend=numHighlights;
    hold on
    ax=zeros(numLegend, 1);

    for iLg=1:numLegend
        ax(iLg)=plot(NaN,NaN,'.', 'MarkerFaceColor', ncolorAll(iLg, :), 'MarkerEdgeColor',  ncolorAll(iLg, :), 'markersize', 20); %plotting invisible points of desired colors
    end
    legendText=strcat('motif#', num2str((1:numLegend).'+offsetHL));
    % legend(ax, legendText, 'location', 'bestoutside', 'NumColumns', 5)
    legend(ax, legendText, 'location', 'bestoutside', 'NumColumns', 1)

    figname=strcat(outputFolderName,'motifsAllOnGraph',num2str(config.is3D), fignamExtnd);

    saveas(gcf,figname)


% figure
% hold on
%     ax=zeros(numLegend, 1);
% 
%     for iLg=1:numLegend
%         ax(iLg)=plot(NaN,NaN,'.', 'MarkerFaceColor', ncolorAll(iLg, :), 'MarkerEdgeColor',  ncolorAll(iLg, :), 'markersize', 20); %plotting invisible points of desired colors
%     end
%     legendText=strcat('motif#', num2str((1:numLegend).'));
%     legend(ax, legendText, 'location', 'bestoutside')

    axis equal
    axis off


    close(gcf);