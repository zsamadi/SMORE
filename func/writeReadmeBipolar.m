function done=writeReadmeBipolar(filename, fixedTypes, gcConfig)



fileID = fopen(filename,'w');




starLine=repelem('*', 1, 50);
starLine=strcat(starLine, '\n');

fprintf(fileID,starLine);
fprintf(fileID,'MATLAB gSTREME\n');
fprintf(fileID,starLine);

fprintf(fileID,'\n');
fprintf(fileID,'motif length: ');
fprintf(fileID,'%d ', gcConfig.W);
fprintf(fileID,'\n');


fprintf(fileID,'\n');
fprintf(fileID,'Fixed cell types: ');
fprintf(fileID,'%d ', fixedTypes);
fprintf(fileID,'\n');

if any(gcConfig.randiHold)
    fprintf(fileID,'\n');
    fprintf(fileID,'positive hold out sections: ');
    fprintf(fileID,'%d ', gcConfig.randiHold);
    fprintf(fileID,'\n');

    fprintf(fileID,'\n');
    fprintf(fileID,'negetive hold out sections: ');
    fprintf(fileID,'%d ', gcConfig.randiHoldR);
    fprintf(fileID,'\n');

else
    fprintf(fileID,'same data for holdout and training \n');
end    

fprintf(fileID,'shuffle mode is %s\n', gcConfig.shuffleMode);   

LogicalStr = {'false', 'true'};

fprintf(fileID,'uniform BG mode is %s\n', LogicalStr{gcConfig.isUBack+1});   
fprintf(fileID,'diff pvalue mode is %s\n', LogicalStr{gcConfig.isDiffMotif+1});   



fprintf(fileID,'\n');

cellTypesOne=gcConfig.cellTypes;


aac=gcConfig.cTypeChars;
aan=string(1:length(aac));

aacn=[aac.',aan.'];

fprintf(fileID,'alphabet letters \n'); 
  fprintf(fileID,'\n');

for iStr=1:length(cellTypesOne)
    strI=aacn(iStr, :);

    fprintf(fileID,'%s\t', strI); 
    fprintf(fileID,'%s\t', cellTypesOne{iStr}); 
    fprintf(fileID,'\n');
end



  fprintf(fileID,'\n');

for iStr=1:length(cellTypesOne)
    strI=aacn(iStr, 1);
    fprintf(fileID,'%s(', strI); 

    fprintf(fileID,'%s:%s)', aacn(iStr, 2),cellTypesOne{iStr}); 
    fprintf(fileID,'\n');
end



fclose(fileID);

done=0;