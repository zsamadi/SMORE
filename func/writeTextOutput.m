function done=writeTextOutput(filename, motifsCell,bkg, textIn,commandText, elapsedTime)

PWM0=bkg.PWM0;


fileID = fopen(filename,'w');
if ~isempty(motifsCell)
    [alength, w]=size(motifsCell{1}.PWM);
else
    alength=1;
    w=1;
end

starLine=repelem('*', 1, 50);
starLine=strcat(starLine, '\n');

fprintf(fileID,starLine);
fprintf(fileID,'MATLAB zSTREME - Sensitive, Thorough, Rapid, Enriched Motif Elicitation\n');
fprintf(fileID,starLine);
fprintf(fileID,'MEME version 5.5.0 (Release date: Wed Sep 7 14:18:26 2022 -0700)\n');
fprintf(fileID,starLine);
fprintf(fileID,'\n');

fprintf(fileID,'\n');
fprintf(fileID,starLine);

fprintf(fileID,'REFERENCE\n');
fprintf(fileID,starLine);

fprintf(fileID,'Timothy L. Bailey,\n');
fprintf(fileID,'"STREME: accurate and versatile sequence motif discovery",\n');

fprintf(fileID,'Bioinformatics, Mar. 24, 2021.\n');
fprintf(fileID,starLine);
fprintf(fileID,'\n');

fprintf(fileID,'\n');
fprintf(fileID,starLine);
fprintf(fileID,'ALPHABET "ABCDEFGHIJKLMNO"\n');
fprintf(fileID,starLine);

fprintf(fileID,'Background letter frequencies\n');

fprintf(fileID,'A %1.5f B %1.5f C %1.5f D %1.5f E %1.5f F %1.5f G %1.5f H %1.5f I %1.5f J %1.5f K %1.5f L %1.5f M %1.5f N %1.5f O %1.5f\n', PWM0);
fprintf(fileID,'\n');

for im=1:length(motifsCell)
    outMotif=motifsCell{im};
    seedi=outMotif.cSeed;
    seediChar=char(seedi+64);

%     secSeedi=outMotif.secSeeds(1, :);
%     secSeediChar=char(secSeedi+64);


    strText=strcat("MOTIF  ", num2str(im), '-', seediChar,' zSTREME-', num2str(im), '\n');
    fprintf(fileID,strText);

    testPvalue=outMotif.testPvalue/log(10);

    ep=floor(testPvalue);
    ei=10^(testPvalue-ep);
    eText=strcat(num2str(ei), 'e', num2str(ep));

    fprintf(fileID,'letter-probability matrix: alength= %d w= %d nsites= %d scrThr=%1.2f , logEvalue=%2.2f, E= %s\n', alength, w, outMotif.nsites,outMotif.scoreThr,outMotif.testPvalue,eText);



        PWM1=outMotif.PWM;


        for ii = 1:size(PWM1,2)
            fprintf(fileID,'%1.7f\t',(PWM1(:,ii)).');
            fprintf(fileID,'\n');
        end
        fprintf(fileID,'\n');

        seedsToPWM=outMotif.seedsToPWM;
        seedsPvalue=outMotif.seedsPvalue;
        seedsHInfo=outMotif.pnsCPV;

        seedsInfo=[seedsToPWM,seedsPvalue];

    fprintf(fileID,'SEEDS involved in the motif, with pvalues\n');

        for ii = 1:size(seedsInfo,1)
            fprintf(fileID,'seed\t');
            fprintf(fileID,'%d\t',(seedsInfo(ii,1:end-1)));
            fprintf(fileID,'counts\t');
            if any(any(seedsHInfo~=0))              
                fprintf(fileID,'%d\t',(seedsHInfo(ii,1:end-1)));
                fprintf(fileID,'train\\hold pvalues\t');
                fprintf(fileID,'%2.2f\t',[seedsInfo(ii,end), seedsHInfo(ii,end)]);
            else
                check=1;
             end

            fprintf(fileID,'\n');
        end

       
fprintf(fileID,'\n');

%     writematrix(PWMSE,filename,'WriteMode','append','Delimiter','tab')

end
fprintf(fileID,starLine);

fprintf(fileID,textIn);
fprintf(fileID,starLine);

commandText=strcat("COMMAND:    ",commandText, '\n');
fprintf(fileID,commandText);
fprintf(fileID,starLine);

pcName = getenv('COMPUTERNAME');

cpuRext=strcat("CPU:		",pcName, '\n');

fprintf(fileID,cpuRext);

fprintf(fileID,starLine);
fprintf(fileID,'FINALTIME:	%5.2f seconds\n', elapsedTime);

fprintf(fileID,starLine);


fclose(fileID);

done=0;