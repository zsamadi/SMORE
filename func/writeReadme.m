function done=writeReadme(filename, gcConfig)



fileID = fopen(filename,'w');




starLine=repelem('*', 1, 50);
starLine=strcat(starLine, '\n');

fprintf(fileID,starLine);
fprintf(fileID,'MATLAB Smore\n');
fprintf(fileID,starLine);

fprintf(fileID,starLine);

commandText=strcat("COMMAND:    ",gcConfig.commandText, '\n');
fprintf(fileID,commandText);
fprintf(fileID,starLine);



fprintf(fileID,'\n');

cellTypesOne=gcConfig.cellTypes;


aac=gcConfig.cTypeChars;
aan=string(1:length(aac));

aacn=[aac(:),aan(:)];

fprintf(fileID,'alphabet letters \n'); 
  fprintf(fileID,'\n');

for iStr=1:length(cellTypesOne)
    strI=aacn(iStr, :);

    fprintf(fileID,'%s\t', strI); 
    fprintf(fileID,'%s\t', cellTypesOne{iStr}); 
    fprintf(fileID,'\n');
end




fclose(fileID);

done=0;