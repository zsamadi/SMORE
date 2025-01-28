function gExHPlot(gExResults, specs)   
[~, sorID]=sort(specs.alphabet);
rankID(sorID)=1:length(sorID);


 % save Results
 % addressC=gExResults.addressC;
 % 
 %     addressAll0=vertcat(addressC{:});
 %     addressObvs=ismember(addressAll0(:, 1), [(1:5).';9]);
 %     gExResults=structfun(@(x) x(~addressObvs), gExResults, 'UniformOutput', false);
 %     addressC=gExResults.addressC;
 % 
 %     addressAll0=vertcat(addressC{:});
 % 
 % 
 %     addressAll0WOP=addressAll0;
 %    addressAll0WOP(:, 2)=[];
 %    [~, aPU]=unique(addressAll0WOP, 'rows');
 % 
 %    aa13=addressAll0(:, 1:3);
 %    [~, ~, aajU]=unique(aa13, 'rows');
 %    aac=accumarray(aajU, 1);
 %    aac=aac(aajU);
 %    aacf=aac>2;
 %    gExResults=structfun(@(x) x(aacf), gExResults, 'UniformOutput', false);




    geneSet=gExResults.geneSet;
    pvalueMedTheoC=gExResults.pvalueMedTheoC;
    foldChangeC=gExResults.foldChangeC;
    addressC=gExResults.addressC;
    medMotC=gExResults.medMotC;
    medtypeC=gExResults.medtypeC;
    numMotCellsC=gExResults.numMotCellsC;



    pvalAllHMapMin=specs.pvalHMapMin;
    dMedAllHMapMin=specs.dMedHMapMin;
    geneAnnotes=specs.geneAnnotes;
    isGEByTissue=specs.isGEByTissue;



    outputFolderName=specs.outputFolderName;


    pvalAllTheoMin=3;
    tablePvalQMin=0;



    %% plot volcano

    addressAll0=vertcat(addressC{:});
    % addressAll0(:, 3)=rankID(addressAll0(:, 3));




    % heatmap for all cases, ordered with address
    lengthV=zeros(length(pvalueMedTheoC),1);

    for iPV=1:length(pvalueMedTheoC)
        lengthV(iPV)=length(pvalueMedTheoC{iPV});

    end



    pvalAllHist=vertcat(pvalueMedTheoC{:});

    pvalAllHist=abs(log(pvalAllHist));

    isInfPval=isinf(pvalAllHist);
    pvalAllHist(isInfPval)=ceil(max(pvalAllHist(~isInfPval))+1);

    % [pvalAllHistSort, iHSort]=sort(pvalAllHist);


    pvalAllTheo=vertcat(pvalueMedTheoC{:});

    pvalAllTheo=abs(log(pvalAllTheo));


    isInfPval=isinf(pvalAllTheo);
    pvalAllTheo(isInfPval)=ceil(max(pvalAllTheo(~isInfPval))+1);



    % pvalAllFitSort=pvalAllTheo(iHSort);


    pvalAllFlag=pvalAllTheo>pvalAllTheoMin;

    geneSetAll=vertcat(geneSet{:});

    geneSetAllS=geneSetAll(pvalAllFlag);
    pvalAllHistS=pvalAllHist(pvalAllFlag);


    pvalAllFitS=pvalAllTheo(pvalAllFlag);

    foldChangeAll=vertcat(foldChangeC{:});




    figure('visible','off')
    scatter(foldChangeAll, pvalAllHist, 'filled')

    grid minor
    legend('motif gene expression pvalues')
    xlabel('delta median')
    ylabel('absolute log pvalue')

    figName='volcano';
    figName=strcat(outputFolderName, figName, '.jpg');

    saveas(gcf,figName);
    %closegcf);



    textingFlag= pvalAllHistS>=quantile(pvalAllFitS, tablePvalQMin);
    foldChangeAllS=foldChangeAll(pvalAllFlag);
    addressAll=repelem(addressAll0, lengthV,1);
    numMotAllV=vertcat(numMotCellsC{:});
    numMotAllV=repelem(numMotAllV, lengthV,1);
    addressAllS=addressAll(pvalAllFlag, :);


    geneTextD=geneSetAllS(textingFlag);
    geneText=geneAnnotes(geneTextD);

    geneName=geneText(:);

    motifAddr=addressAllS(textingFlag, :);
    motifNumber=motifAddr(:, 1);
    motifPosition=motifAddr(:,2);
    cellType=motifAddr(:,3);
    if specs.isGEByTissue
        section_id=motifAddr(:,4);
    end

    pvalue=pvalAllFitS(textingFlag);
    deltaMedian=foldChangeAllS(textingFlag);
    medMotVec=vertcat(medMotC{:});
    medtypeVec=vertcat(medtypeC{:});
    medMotVecS=medMotVec(pvalAllFlag);
    medtypeVecS=medtypeVec(pvalAllFlag);
    motifMedian=medMotVecS(textingFlag);
    typeMedian=medtypeVecS(textingFlag);
    numMotAllVS=numMotAllV(pvalAllFlag, :);
    numCellsMot=numMotAllVS(textingFlag, :);

    numCellMotif=numCellsMot(:, 1);
    numCellType=sum(numCellsMot, 2);






    if specs.isGEByTissue
        asignificantGenesTable=table(geneName,motifNumber,motifPosition,cellType, section_id, pvalue, deltaMedian, motifMedian, typeMedian, numCellMotif, numCellType);
    else
        asignificantGenesTable=table(geneName,motifNumber,motifPosition,cellType, pvalue, deltaMedian, motifMedian, typeMedian, numCellMotif, numCellType);
    end 

    [~, ist]=sort(pvalue, 'descend');

    asignificantGenesTable=asignificantGenesTable(ist, :);



    csvFileName='geaTable';
    csvFileName=strcat(outputFolderName, csvFileName, '.csv');
    if exist(csvFileName, "file")
        delete(csvFileName)
    end

    writetable(asignificantGenesTable,csvFileName)


    %%


    heatMapVal0=zeros(length(pvalueMedTheoC), length(geneAnnotes));
    heatMapValMD0=zeros(length(pvalueMedTheoC), length(geneAnnotes));

    for iPV=1:length(pvalueMedTheoC)
        heatMapValMD0(iPV, geneSet{iPV})=foldChangeC{iPV};

        heatMapVal0(iPV, geneSet{iPV})=abs(log(pvalueMedTheoC{iPV}));
    end




    heatmapYLabel=addressAll0;
    [heatmapYLabel, iAddrSort]=sortrows(heatmapYLabel);
    addressAll0=addressAll0(iAddrSort, :);
    [heatmapYLabelU, iU]=unique(heatmapYLabel(:, 1));
    if isGEByTissue
        heatmapYLabel=strcat(num2str(heatmapYLabel(:, 1)),'-', num2str(heatmapYLabel(:, 2)),'-', num2str(heatmapYLabel(:, 3)),'-', num2str(heatmapYLabel(:, 4)));
    else
        heatmapYLabel=strcat(num2str(heatmapYLabel(:, 1)),'-', num2str(heatmapYLabel(:, 2)),'-', num2str(heatmapYLabel(:, 3)));
    end

    % heatmapYLabel=strcat(num2str(heatmapYLabel(:, 1)),'-', num2str(heatmapYLabel(:, 2)),'-', num2str(heatmapYLabel(:, 3)));

    heatmapYLabelStr=repelem({''},length(heatmapYLabel));

    iUL=iU+[iU(2:end);length(heatmapYLabel)];
    iUL=round(iUL/2);

    for iiU=1:length(iU)

        heatmapYLabelStr{iUL(iiU)}=num2str(heatmapYLabelU(iiU));
    end

    heatMapVal=heatMapVal0(iAddrSort, :);
    heatMapValMD=heatMapValMD0(iAddrSort, :);



    %%%%  heatmap

    f=figure('visible','off');
    f.Position(3:4)=2*f.Position(3:4);

    h=heatmap(heatMapVal,'XLabel' , 'gene','YLabel' , 'motif number','Colormap',cool,'CellLabelColor','none');



    h.GridVisible = 'on';
    c=gray;
    c = flipud(c);
    colormap(c);
    Ax = gca;
    Ax.XDisplayLabels = nan(length(Ax.XDisplayData),1);
    Ax.YDisplayLabels = nan(length(Ax.YDisplayData),1);
    % set(gca,'ColorScaling','log') % Log scale
    % axs = struct(gca);
    % cb = axs.Colorbar;
    %
    % cb.Ticks = unique(floor(linspace(min(heatMapVal(:)), max(heatMapVal(:)), 8))) ;
    %
    %
    % Ax.XDisplayLabels=geneAnnotes;

    % set(struct(h).NodeChildren(3), 'XTickLabelRotation', 0); % put instead of the last example line
    Ax.YDisplayLabels=heatmapYLabel;

    colorbar
    title('log pvalue');
    figName=strcat(outputFolderName, 'allHeatmap.jpg');
    saveas(gcf, figName)
    %closegcf);

    


    %%%%%%%%%%%%%%%%%%


    % heatmap for selected ones


    geneFlags=any(heatMapVal>pvalAllHMapMin & abs(heatMapValMD)>dMedAllHMapMin);
    postFalgs=any(heatMapVal>pvalAllHMapMin  & abs(heatMapValMD)>dMedAllHMapMin, 2);






    % 
    % motifNum=addressAll0(:, 1);
    % numSigs=accumarray(motifNum, postFalgs);
    % postFalgs=numSigs>0;
    % postFalgs=postFalgs(motifNum);
    addressAllSel=addressAll0(postFalgs, :);
    heatMapValSel=heatMapVal(postFalgs, :);
    heatMapValMDSel=heatMapValMD(postFalgs, :);

    heatMapValSel=heatMapValSel(:, geneFlags);
    heatMapValMDSel=heatMapValMDSel(:, geneFlags);

    geneAnnotesSell=geneAnnotes(geneFlags);
    heatmapYLabelSel=heatmapYLabel(postFalgs,:);
    % f=figure('visible','off');
    % f.Position(3:4)=2*f.Position(3:4);
    % 
    % h=heatmap(heatMapValSel,'FontSize',12,'Colormap',cool,'CellLabelColor','none');
    % c=gray;
    % c = flipud(c);
    % colormap(c);
    % % hHeatmap = struct(h).Heatmap;
    % % hHeatmap.GridLineStyle = ':';
    % 
    % Ax = gca;
    % Ax.XDisplayLabels=geneAnnotesSell;
    % Ax.YDisplayLabels=heatmapYLabelSel;
    % set(gca,'ColorScaling','log')
    % axs = struct(gca);
    % cb = axs.Colorbar;
    %
    % cb.Ticks = unique(floor(linspace(min(heatMapVal(:)), max(heatMapVal(:)), 8))) ;


    % figName=strcat(outputFolderName, 'selectedHeatmap.jpg');
    % saveas(gcf, figName)

    % heatmap for selected ones, ordered with type

    heatmapYLabel=addressAllSel;
    if ~isempty(heatmapYLabel)
        sortAddr=[heatmapYLabel(:, end),heatmapYLabel(:, end-2), heatmapYLabel(:, 1:end-2)];
        [~, iAddrSort]=sortrows(sortAddr);
        heatmapYLabel=heatmapYLabel(iAddrSort, :);
        if isGEByTissue
            heatmapYLabel=strcat(num2str(heatmapYLabel(:, 1)),'-', num2str(heatmapYLabel(:, 2)),'-', num2str(heatmapYLabel(:, 3)),'-', num2str(heatmapYLabel(:, 4)));
        else
            heatmapYLabel=strcat(num2str(heatmapYLabel(:, 1)),'-', num2str(heatmapYLabel(:, 2)),'-', num2str(heatmapYLabel(:, 3)));
        end
            

        heatMapValSel=heatMapValSel(iAddrSort, :);
        heatMapValMDSel=heatMapValMDSel(iAddrSort, :);

        [geneAnnotesSell, idGSort]=sort(geneAnnotesSell);
        % [~, idGSort]=sortrows(heatMapValSel.');
        % geneAnnotesSell=geneAnnotesSell(idGSort);

        heatMapValSel=heatMapValSel(:, idGSort);
        heatMapValMDSel=heatMapValMDSel(:, idGSort);

    else
        heatmapYLabel=heatmapYLabelSel;
    end



    f=figure('visible','on');
    f.Position(3:4)=2*f.Position(3:4);
    f.Position(1:2)=1/2*f.Position(1:2);
    % h=heatmap(heatMapValSel,'FontSize',6,'Colormap',cool,'CellLabelColor','white', 'CellLabelFormat', '%.0f');
    h=heatmap(heatMapValSel,'FontSize',6,'Colormap',cool,'CellLabelColor','none');

    c=gray;
    c = flipud(c);
    colormap(c);
    if ~isempty(heatMapValSel)
        clim([0,min(50, max(heatMapValSel(:)))])
    end
    % hHeatmap = struct(h).Heatmap;
    % hHeatmap.GridLineStyle = ':';

    Ax = gca;
    Ax.XDisplayLabels=geneAnnotesSell;
    Ax.YDisplayLabels=heatmapYLabel;
    % set(gca,'ColorScaling','log')
    % axs = struct(gca);
    % cb = axs.Colorbar;
    %
    % cb.Ticks = unique(floor(linspace(min(heatMapVal(:)), max(heatMapVal(:)), 8)));
    title('log pvalue');

    figName=strcat(outputFolderName, 'selectedSortedHeatmap.jpg');
    saveas(gcf, figName)
    %closegcf);

    cgo_all = clustergram(heatMapValSel.*sign(heatMapValMDSel),'Standardize','none', 'linkage', 'complete','RowPDist', 'correlation','ColumnPDist', 'correlation' );
    heatmapYLabelCell=arrayfun(@(i) {heatmapYLabel(i, :)}, (1:size(heatmapYLabel, 1)));
    % set(cgo_all,'RowLabels',heatmapYLabelCell,'ColumnLabels',geneAnnotesSell)
    sortColumn=cgo_all.ColumnLabels;
    sortColumn=sortColumn.';
    sortColumn=cell2mat(sortColumn);
    sortColumn=str2num(sortColumn);

    set(cgo_all,'RowLabels',heatmapYLabelCell,'ColumnLabels',geneAnnotesSell, 'Colormap', redbluecmap)
    f=figure('visible','off');
    % f.Position(3:4)=3*f.Position(3:4);
    % f.Position(1:2)=1/3*f.Position(1:2);

    plot(cgo_all, f)

    figName=strcat(outputFolderName, 'selectedClutergram.jpg');
    saveas(gcf, figName)



    % f=figure('visible','on');
    % f.Position(3:4)=2*f.Position(3:4);
    % f.Position(1:2)=1/2*f.Position(1:2);
    % 
    % % f.Position(3:4)=2*f.Position(3:4);
    % 
    % h=heatmap(heatMapValSel(:, sortColumn),'FontSize',6,'Colormap',cool,'CellLabelColor','none');
    % 
    % c=gray;
    % c = flipud(c);
    % colormap(c);
    % clim([0,min(400, max(heatMapValSel(:)))])
    % 
    % Ax = gca;
    % Ax.XDisplayLabels=geneAnnotesSell(sortColumn);
    % Ax.YDisplayLabels=heatmapYLabel;
    % 
    % 
    % title('log pvalue');
    % 
    % figName=strcat(outputFolderName, 'selectedPvalClHeatmap.jpg');
    % saveas(gcf, figName)



    f=figure('visible','on');
    % f.Position(3:4)=3*f.Position(3:4);
    f.Position(3:4)=2*f.Position(3:4);
    f.Position(1:2)=1/2*f.Position(1:2);

    h=heatmap(heatMapValMDSel,'FontSize',6,'Colormap',cool,'CellLabelColor','white', 'CellLabelFormat', '%1.1f');
    c=colormap(lbmap(256,'BrownBlue'));
    colormap(c);
    cMax=max(abs(heatMapValMDSel(:)));

    cMax=min([cMax, 15]);

    if ~isempty(heatMapValSel)
        clim([-cMax,cMax])
    end


    
    % hHeatmap = struct(h).Heatmap;
    % hHeatmap.GridLineStyle = ':';

    Ax = gca;
    Ax.XDisplayLabels=geneAnnotesSell;
    Ax.YDisplayLabels=heatmapYLabel;

    title('delta median');

    figName=strcat(outputFolderName, 'selectedDMedianHeatmap.jpg');
    saveas(gcf, figName)
    %closegcf);






















