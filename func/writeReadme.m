function done=writeReadme(filename, fixedTypes, gcConfig)



fileID = fopen(filename,'w');




starLine=repelem('*', 1, 50);
starLine=strcat(starLine, '\n');
fprintf(fileID,starLine);
fprintf(fileID,gcConfig.sampleSpecs);
fprintf(fileID,'\n');

fprintf(fileID,starLine);
fprintf(fileID,gcConfig.statURPENG);
fprintf(fileID,'\n');


fprintf(fileID,starLine);
fprintf(fileID,'MATLAB Smore\n');
fprintf(fileID,starLine);

fprintf(fileID,'\n');
fprintf(fileID,'motif length: ');
fprintf(fileID,'%d ', gcConfig.W);
fprintf(fileID,'\n');


fprintf(fileID,'\n');
fprintf(fileID,'Fixed cell types: ');
fprintf(fileID,'%d ', fixedTypes);
fprintf(fileID,'\n');



fprintf(fileID,'\n');

cellTypesOne=gcConfig.cellTypesOne;




aac=gcConfig.cTypeChars;
aan=string(gcConfig.cellTypes);

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