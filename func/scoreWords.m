function [PWMScore, PWMScoreLetter]=scoreWords(seeds, PWMS, scrSpecs)
% Scores the input seeds with PWMS matrix
% Input arguments:
%   - seeds: input seeds to be scored
%   - PWMS   : input PWM score matrix
% Output arguments:
%   - PWMScore: PWM score of the input seeds
%   - PWMScoreLetter: letter by letter score of the input seeds
% used in seqFilterNew

order=scrSpecs.mkvOrder;



[numCells, W]=size(PWMS);

% initialize PWM Score matrix, score for erased letters are set to -10
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
% PWMScore=round(sum(PWMScoreLetter,2), 10);
PWMScore=sum(PWMScoreLetter,2);

PWMScore=gather(PWMScore);
PWMScoreLetter=gather(PWMScoreLetter);


if order>0
    scoreBack=getScoreBack(seeds,scrSpecs);
    
    PWMScoreLetter(~wordsErs, :)=PWMScoreLetter(~wordsErs, :)-scoreBack;
    PWMScore=sum(PWMScoreLetter,2);

end



