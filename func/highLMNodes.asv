function highLMNodes(G, ppHMNodesCell,outputFolderName, config)
fignamExtnd=config.fignamExtnd;

cTypeChars=config.cTypeChars;
cPTypes=(G.Nodes.label(:, 1)).';

xcoordsTotal=G.Nodes.Coordinates(:, 1);
ycoordsTotal=G.Nodes.Coordinates(:, 2);

if config.is3D
    zcoordsTotal=G.Nodes.Coordinates(:, 3);
end

numExtMotifs=length(ppHMNodesCell);
W=config.W;

for iMG=1:numExtMotifs

    figure('visible','off');

    if config.is3D
        h3=plot(G,'XData',xcoordsTotal,'YData',ycoordsTotal, 'ZData',zcoordsTotal);
    else
        h3=plot(G,'XData',xcoordsTotal,'YData',ycoordsTotal);
    end



    motifPaths=ppHMNodesCell{iMG};
    cMTypes=cPTypes(motifPaths);

    [cMTypesU, iu, ju]=unique(cMTypes, 'rows');
    cCounts=accumarray(ju, 1);
    [cCounts, ist]=sort(cCounts, 'descend');
    cMTypesU=cMTypesU(ist, :);
    cCountsTotal=sum(cCounts);
    ncolor=[1 0 1;0 1 0;1 0.64 0;1 0 0;0 0 0;0 1 1;rand(length(iu),3)];



    for iCTU=1:length(iu)
        cli=ncolor(iCTU, :);
        seedPaths=motifPaths(ju==ist(iCTU), :);
        seedPaths=seedPaths(:);


        highlight(h3,seedPaths,'NodeColor',cli, 'EdgeColor',cli);
    end


    titleStr=sprintf("motif %d- total Counts %d", iMG, cCountsTotal);
    for iStr=1:min(length(cCounts), W+1)
        countStr=sprintf(":%d ", cCounts(iStr));
        titleStr=strcat(titleStr, countStr);
    end
    title(titleStr)

    grid minor
    % 
    % 
    % 
    % set(gca,'xtick',(sectionU-1)*15000,'xticklabel','')
    % set(gca,'ytick',(AnimalIDIdx-1)*10000,'yticklabel','')


    xlabel('X')
    ylabel('Y')
    numLegend=min(6, length(iu));
    hold on
    ax=zeros(numLegend, 1);

    for iLg=1:numLegend
        ax(iLg)=plot(NaN,NaN,'.', 'MarkerFaceColor', ncolor(iLg, :), 'MarkerEdgeColor',  ncolor(iLg, :), 'markersize', 20); %plotting invisible points of desired colors
    end
    legendText=cTypeChars(cMTypesU(1:numLegend, :));
    legend(ax, legendText)



    figname=strcat(outputFolderName,'motifGroupOnGraph', num2str(iMG),num2str(W),num2str(config.is3D), fignamExtnd);

    saveas(gcf,figname)
    close(gcf);

end