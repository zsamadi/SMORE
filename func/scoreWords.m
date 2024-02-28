function [PWMScore, PWMScoreLetter]=scoreWords(seeds, PWMS, scrSpecs)


order=scrSpecs.mkvOrder;


[numCells, W]=size(PWMS);

PWMScoreLetter=-10*ones(size(seeds));

wordsErs=any(seeds==0, 2);
if any(wordsErs)
    seeds=seeds(~wordsErs, :);
end

seeds=gpuArray(seeds);
PWMS=gpuArray(PWMS);


wordsC=seeds+(0:W-1)*numCells;
wordsC=wordsC(:);
PWMScore=PWMS(wordsC);
PWMScoreLT=reshape(PWMScore, size(seeds,1), W);
if any(wordsErs)
    PWMScoreLetter(~wordsErs, :)=PWMScoreLT;
else
    PWMScoreLetter=PWMScoreLT;
end
PWMScore=sum(PWMScoreLetter,2);

PWMScore=gather(PWMScore);
PWMScoreLetter=gather(PWMScoreLetter);


if order>0
    scoreBack=getScoreBack(seeds,scrSpecs);

    PWMScoreLetter(~wordsErs, :)=PWMScoreLetter(~wordsErs, :)-scoreBack;
    PWMScore=sum(PWMScoreLetter,2);

end



