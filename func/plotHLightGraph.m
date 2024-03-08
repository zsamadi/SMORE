function h=plotHLightGraph(G, cellSubtypeVec,options)

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
[~, ist]=sort(cellCnts, 'descend');

cellSubtypeVecU=cellSubtypeVecU(ist);
cellSubtypeVecU=cellSubtypeVecU(:);



    numHighlight=min(10, length(cellSubtypeVecU));
    figure('visible','off');
    ncoloro=[0 0.4470 0.7410;0.9290 0.6940 0.1250;0.8500 0.3250 0.0980;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840;0 0 0;0 1 1;1 0 1;0 1 0;1 0.64 0;1 0 0;rand(length(cellSubtypeVecU),3)];


    plotEdges=options.isPlotEdges;

    edgeClrStr=["none", "b"];
    edgeClr=edgeClrStr(double(plotEdges)+1);


    if options.is3D


                ax1=subplot(4, 3, [1,2,4,5,7,8]); 

                h=plot(ax1,G,'XData',xcoordsTotal,'YData',ycoordsTotal,'ZData',zcoordsTotal,'NodeColor',[0.8    0.8    0.8], 'EdgeColor',edgeClr);

    else
            h=plot(G,'XData',xcoordsTotal,'YData',ycoordsTotal, 'EdgeColor',edgeClr, 'MarkerSize',1);
    end

        for iu=1:length(cellSubtypeVecU)

        if iu<=numHighlight
            highlight(h,nodesAll(cellSubtypeVec==cellSubtypeVecU(iu)),'NodeColor',ncoloro(iu, :))
        else
            highlight(h,nodesAll(cellSubtypeVec==cellSubtypeVecU(iu)),'NodeColor',[0.6980    0.7451    0.7098])
        end


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
            legendText=strcat(options.ctAnnot(cellSubtypeVecU(1:numHighlight)), '(', string(cellSubtypeVecU(1:numHighlight)), ')');
    else

             legendText=string(cellSubtypeVecU(1:numHighlight));
    end
    
    if numHighlight<length(cellSubtypeVecU)
        ax(end)=plot(NaN,NaN,'.', 'MarkerFaceColor', [0.6980    0.7451    0.7098], 'MarkerEdgeColor',[0.6980    0.7451    0.7098], 'markersize', 20); %plotting invisible points of desired colors
        legendText=[legendText;"others"];

    end


    numLegColumns=ceil(numHighlight/20);


    legend(ax,legendText,'NumColumns',numLegColumns, 'Location', 'bestoutside')




    if options.isShuffled>0
         figname=strcat(folderName, 'GraphWithShuffleTypes.jpg');
         titleTex='shuffled';

    else

        figname=strcat(folderName, 'GraphWithTypes.jpg');
        titleTex='primary';


    end
    title(titleTex)

    hold off
    if options.is3D
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

 if ~options.is3D
    saveas(gcf,figname)
 end


   