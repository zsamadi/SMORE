
function [xLDimOut, geneExpression, cellSubtype] =getPCA(cellTable, options)
    offset=options.offset;
    cellSubtype=cellTable.cellSubtypeOrig;

    
    % gene expresion for the naive cells of the selected nodes
    geneExpression=cellTable(:, 10:end-offset-21);
    geneExpression=table2array(geneExpression);
    
    % putting aside genes with nan expresion
    nanExpresionFlag = isnan(geneExpression);
    nanExpresionFlagColumn=any(nanExpresionFlag, 1);
    geneExpression(:, nanExpresionFlagColumn)=[];

    aa=cellTable.Properties.VariableNames;
    iidl=[];
    for ii=1:length(aa)
        if strcmpi(aa{ii}(1:end-1), 'blank_')
            iidl=[iidl, ii];
        end
    end
    geneExpression(:,iidl)=[];




    % geneExpression=cellTable(:, end-offset-21: end-offset);
    % geneExpression=table2array(geneExpression);
    % 
    % nanExpresionFlag = isnan(geneExpression);
    % nanExpresionFlagColumn=any(nanExpresionFlag, 1);
    % geneExpression(:, nanExpresionFlagColumn)=[];
    % 
    % 
    % zeroCells=sum(geneExpression, 2)==0;
    % 
    % geneExpression(zeroCells, :)=[];

    % geneExpressionN0=geneExpression(geneExpression>0);
    % geneExpressionN0Min=min(geneExpressionN0);
    % geneExpression=10*geneExpression/geneExpressionN0Min;


    geneExpression21=cellTable(:, end-offset-21+1: end-offset);
    geneExpression21=table2array(geneExpression21);

    nanExpresionFlag = isnan(geneExpression21);
    nanExpresionFlagColumn=any(nanExpresionFlag, 1);
    geneExpression21(:, nanExpresionFlagColumn)=[];

    geneExpression21N0=geneExpression21(geneExpression21>0);
    geneExpression21N0Min=min(geneExpression21N0);
    % geneExpression21=10*geneExpression21/geneExpression21N0Min;

    geneExpression=[geneExpression,geneExpression21];

    zeroCells=sum(geneExpression, 2)==0;

    geneExpression(zeroCells, :)=[];

    cellSubtype(zeroCells)=[];


    
    % geneExpression=10000*geneExpression./sum(geneExpression, 2);
    
    %%% Not A Clusters
    
    % to see if clustering was correct, check it on non A cell types
    
    % typExcld=[];
    geneExUMAP=geneExpression;
    geneExUMAPSum=sum(geneExUMAP, 2);
    geneExUMAPSum=geneExUMAPSum/median(geneExUMAPSum);

    % geneExUMAP=geneExUMAP./geneExUMAPSum;

    geneExUMAP=log(1+geneExUMAP);

    % geneExpression = (mapstd(geneExpression.')).';  % Normalize data




    
    % geneExpression=(geneExpression-mean(geneExpression))./std(geneExpression);
    % geneExUMAP=(geneExUMAP-mean(geneExUMAP))./std(geneExUMAP);

if options.isDoPCA   
    x = mapstd(geneExUMAP');  % Normalize data
    xLDim = processpca(x); 
    % [coeff,score,latent,tsquared,explained,mu] = pca(x.');


    xLDim=xLDim.';

    if strcmpi(options.demReduct, 'umap')    
        args= {'verbose' ,'none', 'fast_approximation' , true};    
        xLDim=run_umap(xLDim(:, 1:20), 'cluster_output', 'none', 'cluster_detail', 'medium',  args{:});
    else
        xLDim=xLDim(:, 1:options.pcaDim);

    end

    xLDimOut=zeros(length(zeroCells), size(xLDim, 2));
    xLDimOut(~zeroCells, :)=xLDim;
else
    xLDimOut=[];
end
    check=1;

