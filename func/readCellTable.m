function [T, cellTypesOut, ambgTypes]=readCellTable(folderName, options)


% csvName='Hypothalamus';
% filename=strcat(folderName,'\',csvName, '.csv');
% T = readtable(filename);
filename=strcat(folderName, 'THP.mat');

load(filename, 'T')

C1=T.Neuron_cluster_ID;
[C1U, ~, jU]=unique(C1);
if options.isJClusterID
        C1ULen=length(C1U);
        c1UN=cell(C1ULen, 1);
    for ii=1:C1ULen
    
        tmp=C1U{ii};
        tmpL=length(tmp);
        cL=min(tmpL, 3);
        c1UN{ii}=tmp(1:cL);
    end
        
        T.Neuron_cluster_ID=c1UN(jU);
        
        [C1U, ~, jU]=unique(T.Neuron_cluster_ID);
else
    C1ULen=length(C1U);
    
    for ii=1:C1ULen
        tmpL=length(C1U{ii});
        % if strcmpi(C1U{ii}, 'E-14')
        %     tmpL=tmpL-1;
        % elseif strcmpi(C1U{ii}, 'E-1')
        %     tmpL=tmpL+1;
        % end




        C1U{ii}=strcat(num2str(tmpL), C1U{ii});
    end
    
    Neuron_cluster_IDT=C1U(jU);
    [C1U, ~, jU]=unique(Neuron_cluster_IDT);

end






C1UDg=(1:length(C1U)).';
C1UDg=C1UDg(jU);


C2=T.Cell_class;
[C2U, ~, jU]=unique(C2);
C2UDg=(1:length(C2U)).';
C2UDg=C2UDg(jU);
% 
% if options.isJClusterID
%     C12Dg=C1UDg;
% else
%     C12Dg=[C1UDg, C2UDg];
% end  


C12Dg=[C1UDg, C2UDg];

[cellSubClusterVecU, iU, jU]=unique(C12Dg, 'rows');
cellSubClusterVecU=(1:length(cellSubClusterVecU)).';

cellSubtypeVec=cellSubClusterVecU(jU);
T.cellSubtypeOrig=cellSubtypeVec;

% cSTCnt=accumarray(jU, 1);
% 
% [~, idx]=sort(cSTCnt, 'descend');
% cellSubClusterVecU=cellSubClusterVecU(idx);
% 
% idxAll=(1:length(cellSubtypeVec));
% 
% if options.ambigRatio>0  
%     % if options.isRandAmbig
%     %     % rng(1);
%     %     randPerm20=randperm(20);
%     %     cellSubClusterVecU= cellSubClusterVecU([randPerm20, (21:length(cellSubClusterVecU))]);
%     % end
% 
%     cellSubClusterVecAmbig=options.ambgSWTs;
% 
%     if options.ambigType==0
%         ambigType0=cellSubClusterVecU(end);
%         ambigType1=cellSubClusterVecU(end-1);
%         cellSubClusterVecU1=cellSubClusterVecU;
%         cellSubClusterVecU1(ismember(cellSubClusterVecU1, cellSubClusterVecAmbig))=[];
%         for iST=1:length(cellSubClusterVecU1(1:end-2))
%             ctFlag=(cellSubtypeVec==cellSubClusterVecU1(iST));
%             idxi=idxAll(ctFlag);
%             idxPermi=randperm(length(idxi));
%             idxi=idxi(idxPermi);
%             switchNumberi=ceil(options.ambigRatio1*length(idxi));
%             cellSubtypeVec(idxi(1:switchNumberi))=ambigType1;
%         end
% 
%         ambgTypes=[ambigType0, ambigType1];
%     else
%         ambigType0=options.ambigType;
%         ambgTypes=ambigType0;
% 
%     end
% 
% 
% 
%     for iST=1:length(cellSubClusterVecAmbig)
%         ctFlag=(cellSubtypeVec==cellSubClusterVecAmbig(iST));
%         idxi=idxAll(ctFlag);
%         idxPermi=randperm(length(idxi));
%         idxi=idxi(idxPermi);
%         switchNumberi=ceil(options.ambigRatio*length(idxi));
%         cellSubtypeVec(idxi(1:switchNumberi))=ambigType0;
%     end
% end
% 





% cellNumber=length(cellSubtypeVec);
% 
% cellIndexRand=randperm(cellNumber);
% 
% switchNumber=ceil(options.ambigRatio*cellNumber);
% 
% cellSubtypeVec(cellIndexRand(1:switchNumber))=options.ambigType;


T.cellSubtype=cellSubtypeVec;


cellTypesOut=[T.Neuron_cluster_ID, T.Cell_class];

cellTypesOut=cellTypesOut(iU, :);

if all(options.bregID>0)
    BregmaU=unique(T.Bregma);
    T=T(ismember(T.Bregma, BregmaU(options.bregID)), :);
else

    if options.isNBreg && ~options.isPBreg
        T=T(T.Bregma<0, :);
    elseif options.isPBreg && ~options.isNBreg
        T=T(T.Bregma>0, :);
    end

end




if options.isFemale && ~options.isMale
    T=T(cellfun(@(x) x(1),T.Animal_sex)=='F', :);
elseif options.isMale && ~options.isFemale
    T=T(cellfun(@(x) x(1),T.Animal_sex)=='M', :);
end


if options.isNaive && ~options.isStimu
    T=T(cellfun(@(x) x(1),T.Behavior)=='N', :);
elseif options.isStimu && ~options.isNaive
    T=T(cellfun(@(x) x(1),T.Behavior)~='N', :);
end



cellSubtypeVec=T.cellSubtype;


[cellSubClusterVecU, ~, jU]=unique([cellSubtypeVec;cellSubClusterVecU]);



cSTCnt=accumarray(jU, 1);

[~, idx]=sort(cSTCnt, 'descend');
cellSubClusterVecU=cellSubClusterVecU(idx);

idxAll=(1:length(cellSubtypeVec));
ambgTypes=0;
if options.ambigRatio>0  
    % if options.isRandAmbig
    %     % rng(1);
    %     randPerm20=randperm(20);
    %     cellSubClusterVecU= cellSubClusterVecU([randPerm20, (21:length(cellSubClusterVecU))]);
    % end

    cellSubClusterVecAmbig=options.ambgSWTs;

    ambigType=cellSubClusterVecU(end);




    switchNumberT=0;
    for iST=1:length(cellSubClusterVecAmbig)
        ctFlag=(cellSubtypeVec==cellSubClusterVecAmbig(iST));
        idxi=idxAll(ctFlag);
        idxPermi=randperm(length(idxi));
        idxi=idxi(idxPermi);
        switchNumberi=ceil(options.ambigRatio*length(idxi));
        cellSubtypeVec(idxi(1:switchNumberi))=ambigType;
        switchNumberT=switchNumberi+switchNumberT;
    end
    ambgTypes=ambigType;

    if options.swRatio>0
       swType=cellSubClusterVecU(end-1);

        % cellSubClusterVecU1=cellSubClusterVecU;
        % cellSubClusterVecU1(ismember(cellSubClusterVecU1, cellSubClusterVecAmbig))=[];
        idxAll(ismember(cellSubtypeVec, cellSubClusterVecAmbig))=[];
        switchNumberT=options.swRatio*length(idxAll);
        idxi=randsample(idxAll, switchNumberT);
        cellSubtypeVec(idxi)=swType;

        ambgTypes=[ambigType, swType];
    end

end

T.cellSubtype=cellSubtypeVec;

