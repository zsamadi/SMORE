function h=plotHLightGraph(G, cellSubtypeVec,options)

alphabet=options.alphabet(:);

folderName=options.folderName;

xcoordsTotal=G.Nodes.Coordinates(:, 1);
ycoordsTotal=G.Nodes.Coordinates(:, 2);



if options.is3D
    zcoordsTotal=G.Nodes.Coordinates(:, 3);
end

nodesAll=(1:length(cellSubtypeVec));



rng(1)

[cellSubtypeVecU, ~, jU]=unique(cellSubtypeVec);

cellCnts=accumarray(jU, 1);
[alphabetS, ist]=sort(alphabet(cellSubtypeVecU));

cellSubtypeVecU=cellSubtypeVecU(ist);
cellSubtypeVecU=cellSubtypeVecU(:);



% for embryo3
% cellSubtypeVecU(4)=[];
%
% cellSubtypeVecU(9)=14;
% cellSubtypeVecU(10)=10;
%
% cellSubtypeVecU(11)=21;


numHighlight=min(10, length(cellSubtypeVecU));
% nodesColor=nodesAll(ismember(cellSubtypeVec,cellSubtypeVecU(1:numHighlight)));
% G=subgraph(G, nodesColor);
% 
% cellSubtypeVec=G.Nodes.label(:, 1);
% 
% xcoordsTotal=xcoordsTotal(nodesColor);
% 
% ycoordsTotal=ycoordsTotal(nodesColor);



% figure('visible','off');

ncoloro=[0.9290 0.6940 0.1250;0.8500 0.3250 0.0980;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;0 0 0;0 1 1;1 0 1;0 1 0;1 0.64 0;1 0 0;rand(length(cellSubtypeVecU),3)];



ncoloro= {
    'FF006D';
    'FF870A';
    'FFCE55';
    'FFA900';
    'FF6E00';
    'FF0E48';
    'FFB3DC';
    'F96CBB';
    'AF9DD4';
    '00D9E4';
    '00D9F4';
    '00ACE4';
    'DA84C6';
    '4F98D5';
    '00CBC7'};
ncoloroP=hex2rgb(ncoloro);

% ncoloro=ncoloro(end:-1:1, :);

% 0.8, 1,0;


ncoloro=[0.8500,0.3250,0.0980;0.3, 0.74, 0.93;0.8,1,0;0.46, 0.67, 0.18;0,0.44,0.74;0.63,0.07,0.18;1,0,1;0.49,0.18,0.55;0.92,0.69,0.12];

shuffleColor=[1     7      5     2     4    8     6     9     3];

ncoloro=ncoloro(shuffleColor, :);
ncoloro=[ncoloro;ncoloroP];
ncolorL=options.ncolorL;
% ncolorL= {'3D409A';
% '8183BF';
% 'AC5BA4';
% 'FFC907';
% 'F7EC18';
% 'BD9E6D';
% '69BD45';
% 'ED2024';
% '6ECCDD'};
ncolorL=hex2rgb(ncolorL);
ncoloro=[ncolorL;ncoloro];
% for inl=1:1
%     ncoloro(2, :)=[];
% end

% ncoloro=ncoloro(randperm(size(ncoloro, 1)), :);


plotEdges=options.isPlotEdges;

% edgeClrStr=["none", "b"];
% edgeClr=edgeClrStr(double(plotEdges)+1);
grColor=[0.8,0.8,0.8];
offOn={'off', 'on'};
f=figure('Visible',offOn{options.iShowPlot+1});
    

if options.is3D

    if options.iSideView
        ax1=subplot(4, 3, [1,2,4,5,7,8]);
    else
        f.Position(3:4)=f.Position(3:4)*1.5;
        f.Position(1:2)=f.Position(1:2)*0.75;
        ax1=gca;
    end

    if plotEdges
        h=plot(ax1,G,'XData',xcoordsTotal,'YData',ycoordsTotal,'ZData',zcoordsTotal,'NodeColor',grColor, 'EdgeColor',grColor);
    else
        h=plot(ax1,G,'XData',xcoordsTotal,'YData',ycoordsTotal,'ZData',zcoordsTotal,'NodeColor',grColor, 'EdgeColor','none', 'MarkerSize',2);
    end
else
    f.Position(3:4)=f.Position(3:4)*1.5;
    f.Position(1:2)=f.Position(1:2)*0.75;

    if plotEdges
        h=plot(G,'XData',xcoordsTotal,'YData',ycoordsTotal,'NodeColor',grColor, 'EdgeColor',grColor, 'MarkerSize',1);
    else
        h=plot(G,'XData',xcoordsTotal,'YData',ycoordsTotal,'NodeColor',grColor, 'EdgeColor','none', 'MarkerSize',1);
    end


end

for iu=1:numHighlight

    highlight(h,nodesAll(cellSubtypeVec==cellSubtypeVecU(iu)),'NodeColor',ncoloro(iu, :))



    % if iu<=numHighlight
    %     highlight(h,nodesAll(cellSubtypeVec==cellSubtypeVecU(iu)),'NodeColor',ncoloro(iu, :))
    % else
    %     highlight(h,nodesAll(cellSubtypeVec==cellSubtypeVecU(iu)),'NodeColor',[0.6980    0.7451    0.7098])
    % end


end

hold on

if numHighlight<length(cellSubtypeVecU)
    ax=zeros(numHighlight+1, 1);
else
    ax=zeros(numHighlight, 1);
end

for iu=1:numHighlight
    ax(iu)=plot(NaN,NaN,'.', 'MarkerFaceColor', ncoloro(iu, :), 'MarkerEdgeColor', ncoloro(iu, :), 'markersize', 20); %plotting invisible points of desired colors
end

if ~isempty(options.ctAnnot)
    % legendText=strcat(options.ctAnnot(cellSubtypeVecU(1:numHighlight)), '(', string(cellSubtypeVecU(1:numHighlight)), ')');
    legendText=strcat(options.ctAnnot(cellSubtypeVecU(1:numHighlight)), '(', alphabet(cellSubtypeVecU(1:numHighlight)), ')');

else

    legendText=string(cellSubtypeVecU(1:numHighlight));
end

if numHighlight<length(cellSubtypeVecU)
    ax(end)=plot(NaN,NaN,'.', 'MarkerFaceColor', grColor, 'MarkerEdgeColor',grColor, 'markersize', 20); %plotting invisible points of desired colors
    legendText=[legendText;"others"];

end


numLegColumns=ceil(numHighlight/20);


legend(ax,legendText,'NumColumns',numLegColumns, 'Location', 'bestoutside')




if options.isShuffled>0
    figname=strcat(folderName, 'GraphWithShuffleTypes.jpg');
    titleTex='shuffled tissue map';

else

    figname=strcat(folderName, 'GraphWithTypes.jpg');
    titleTex='tissue map';

end
if iscell(options.tissueIDs)
    for it=1:length(options.tissueIDs)
        titleTex=strcat(titleTex, ',',options.tissueIDs{it});
    end
end


title(titleTex)

hold off
if options.is3D && options.iSideView
    ax2=subplot(4, 3, 3);
    ax3 = subplot(4, 3, 6);
    ax4 = subplot(4, 3, 9);
    ax5 = subplot(4, 3, 10);
    ax6 = subplot(4, 3, 11);
    ax7 = subplot(4, 3, 12);


    xlabel(ax1, 'x');
    ylabel(ax1, 'y');
    zlabel(ax1, 'z');

    copyobj(h,ax2);
    copyobj(h,ax3);
    copyobj(h,ax4);
    copyobj(h,ax5);
    copyobj(h,ax6);
    copyobj(h,ax7);

    xlabel(ax2, 'x');
    ylabel(ax2, 'y');
    zlabel(ax2, 'z');

    xlabel(ax3, 'x');
    ylabel(ax3, 'y');
    zlabel(ax3, 'z');

    xlabel(ax4, 'x');
    ylabel(ax4, 'y');
    zlabel(ax4, 'z');

    xlabel(ax5, 'x');
    ylabel(ax5, 'y');
    zlabel(ax5, 'z');

    xlabel(ax6, 'x');
    ylabel(ax6, 'y');
    zlabel(ax6, 'z');


    xlabel(ax7, 'x');
    ylabel(ax7, 'y');
    zlabel(ax7, 'z');


    view(ax2,90, 0); title(ax2,'view along -x');

    view(ax3,270, 0); title(ax3,'view along +x');

    view(ax4,0, 0); title(ax4,'view along +y');

    view(ax5,0, 270); title(ax5,'view along +z');

    set(ax5, 'Ydir', 'reverse')

    view(ax6,0, 90); title(ax6,'view along -z');

    view(ax7,0, 180); title(ax7,'view along -y');
    set(ax7, 'Zdir', 'reverse')
end

axis equal
% saveas(gcf,figname)

check=1;

if ~options.is3D
    saveas(gcf,figname)
elseif ~options.iSideView
    saveas(gcf,figname)
end

% figure
%
% hold on
%
% if numHighlight<length(cellSubtypeVecU)
%     ax=zeros(numHighlight+1, 1);
% else
%     ax=zeros(numHighlight, 1);
% end
%
% for iu=1:numHighlight
%     ax(iu)=plot(NaN,NaN,'.', 'MarkerFaceColor', ncoloro(iu, :), 'MarkerEdgeColor', ncoloro(iu, :), 'markersize', 20); %plotting invisible points of desired colors
% end
%
% if ~isempty(options.ctAnnot)
%     legendText=strcat(options.ctAnnot(cellSubtypeVecU(1:numHighlight)), '(', (options.alphabet(cellSubtypeVecU(1:numHighlight))).', ')');
% else
%     legendText=string(cellSubtypeVecU(1:numHighlight));
% end
% if numHighlight<length(cellSubtypeVecU)
%     ax(end)=plot(NaN,NaN,'.', 'MarkerFaceColor', grColor, 'MarkerEdgeColor',grColor, 'markersize', 20); %plotting invisible points of desired colors
%     legendText=[legendText;"others"];
%
% end
%
% numLegColumns=ceil(numHighlight/20);
%
%
% legend(ax,legendText,'NumColumns',1, 'Location', 'bestoutside')


