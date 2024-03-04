function [W, hfig] = seqlogoGen(motif,cTypeChars, varargin)
%SEQLOGO displays sequence logos for DNA and protein sequences
%
%   SEQLOGO(SEQS) displays the sequence logo for a set of aligned
%   sequences, SEQS. The logo graphically displays the sequence
%   conservation at a particular position in an alignment of sequences,
%   SEQS, as measured in bits. The maximum sequence conservation per site
%   is log2(4) bits for DNA/RNA and log2(numAlphabet) bits for proteins. SEQS is a
%   string vector or cell array of character vectors containing aligned
%   sequences.
%
%   W = SEQLOGO(SEQS) returns a cell array of a unique symbol list in
%   SEQS, and the information weight matrix used for graphically
%   displaying the logo.
%
%   W = SEQLOGO(...,'DISPLAYLOGO',TF) displays the sequence logo of
%   SEQS when TF is TRUE. The default is TRUE.
%
%   SEQLOGO(...,'ALPHABET',A) specifies that SEQS consists of nucleotides
%   ('NT') or amino acids ('AA'). The default is 'NT'.
%
%   SEQLOGO(...,'STARTAT',STARTPOSITION) specifies the starting position
%   for the sites of interest in SEQS. The default starting position is 1.
%
%   SEQLOGO(...,'ENDAT',ENDPOSITION) specifies the position for the
%   sites of interest in SEQS. The default ending position is the maximum
%   length in SEQS.
% 
%   SEQLOGO(...,'SSCORRECTION',false) specifies not to apply the small
%   sample correction. When there are only a few sample sequences a
%   straightforward calculation tends to overestimate the conservation. By
%   default, SEQLOGO compensates by applying an approximate correction
%   based on the number of sequences. This correction is negligible when
%   more than 50 sequences are used.
% 
%   SEQLOGO(P) displays a sequence logo for P, a sequence profile generated
%   by SEQPROFILE. P is a matrix of size [numAlphabet x seq Length] with the
%   frequency distribution of amino acids. For nucleotides, P is of size 
%   [4 x seq Length] and the DNA alphabet is assumed. P may have 21 (or 5)
%   rows if gaps were included, but SEQLOGO ignores gaps. The sequence
%   conservation is computed without small sample correction. When P
%   contains weighted profiles or symbols counts the profile columns are
%   normalized to the maximum column sum of the profile.
% 
%   [W, H] = SEQLOGO(...) returns the handle of the figure.
% 
%   Example:
% 
%       S = {'ATTATAGCAAACTA',...
%            'AACATGCCAAAGTA',...
%            'ATCATGCAAAAGGA'}
%       % Display the sequence logo of S
%       seqlogo(S)
%
%       % Note that the small sample correction prevents you from seeing
%       % columns with information equal to log2(4) = 2 bits, however you
%       % can also turn this adjustment off:
%       seqlogo(S,'sscorrection',false)
% 
%       % Amino acid sequences
%       S1 = {'LSGGQRQRVAIARALAL'; 
%             'LSGGEKQRVAIARALMN'; 
%             'LSGGQIQRVLLARALAA';
%             'LSGGERRRLEIACVLAL'; 
%             'FSGGEKKKNELWQMLAL'; 
%             'LSGGERRRLEIACVLAL'};
%       seqlogo(S1, 'alphabet', 'aa', 'startAt', 2, 'endAt', 10)
% 
%   Reference:
% 
%       Schneider, T.D., Stephens, R.M., "Sequence Logos: A New Way to
%       Display Consensus Sequences," Nucleic Acids Research, 18, pp.
%       6097-6100, 1990. 
% 
%   See also SEQCONSENSUS, SEQDISP, SEQPROFILE.

%   SEQLOGO(...,'TOFILE',FILENAME) saves the sequence logo in PNG format
%   to a file named FILENAME.png. If FILENAME has no extension, .png is
%   assumed.

%   Copyright 2003-2012 The MathWorks, Inc.

% Validate input data

P=motif.PWM;
background=motif.background;
cSeedStr = strjoin(string(motif.cSeed), ':');

PWMBg=background{1};
PWMBg=PWMBg(:);
PWMBg=ones(size(PWMBg));
PWMBg=PWMBg/sum(PWMBg);

% maxEntropy=-sum(PWMBg.*log2(PWMBg));
maxEntropy=4;


if nargin > 0
    P = convertStringsToChars(P);
end

if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

bioinfochecknargin(nargin, 1, mfilename)

startPos = 1;
endPos = 1;
isAA = true;
NTFlag = false;
displayLogo = true;

numAlphabet=size(motif.PWM, 1);

nSymbols = numAlphabet; % number of symbols in alphabet for NT, numAlphabet for AA 
fileExt = 'png';     % Image file extension
corrError = true;

% cTypeChars = ['ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz', char((198:229))];

cTypeChars=cTypeChars(1:nSymbols);


% cTypeChars = ['ABCDEFGHIJKLMNabcdefhijklmnorstuvwxyzpqgOPQRSTUVWXYZ', char((198:229))];

% 
[~, iSt]=sort(background{1}, 'descend');
% 
% cTypeChars(iSt)=cTypeChars;

ntSymbols = 'ACGT';
symbolList = char(ntSymbols(:));

if nargin > 2
    if rem(nargin,2)== 0
        error(message('bioinfo:seqlogo:IncorrectNumberOfArguments', mfilename));
    end
    okargs = {'alphabet','startat','endat','displaylogo','tofile', 'sscorrection'};
    for j=1:2:nargin-2
        [k, pval] = bioinfoprivate.pvpair(varargin{j}, varargin{j+1}, okargs, mfilename);
        switch(k)
            case 1  % alphabet
                if strcmpi(pval,'aa') % If sequences are Amino Acid
                    isAA = true;
                    nSymbols = numAlphabet; % numAlphabet amino acids
                    symbolList = cTypeChars(:);
                elseif ~strcmpi(pval,'nt') % If sequences are nucleotides
                    warning(message('bioinfo:seqlogo:UnknownAlphabet', upper( pval )));
                    NTFlag = true;
                else
                    NTFlag = true;
                end
            case 2  % start position
                if ~isnumeric(pval) || ~isscalar(pval)
                    error(message('bioinfo:seqlogo:StartAtNotSingleNumericValue'));
                elseif (pval <= 0)
                    error(message('bioinfo:seqlogo:invalidStartAtValue'));
                else
                    startPos = pval;
                end
                
            case 3  % end position
                if ~isnumeric(pval) || ~isscalar(pval)
                    error(message('bioinfo:seqlogo:EndAtNotSingleNumericValue'));
                elseif (pval <= 0)
                    error(message('bioinfo:seqlogo:invalidEndAtValue'));
                else
                    endPos = pval;
                end
            case 4  % showlogo
                displayLogo = bioinfoprivate.opttf(pval);
                if isempty(displayLogo)
                    error(message('bioinfo:seqlogo:InputOptionNotLogical', upper( char( okargs( k ) ) )));
                end
            case 5 % to file
                fileName = pval;
                idx = findstr(fileName,'.');
                if( idx > 0)
                    if ~strcmpi(fileName(idx+1:end), fileExt)
                        warning(message('bioinfo:seqlogo:WrongFileExtension', upper( fileExt ), fileExt));
                        fileName = [fileName(1:idx), fileExt];
                    end
                else
                    fileName = [fileName, '.',  fileExt]; %#ok
                end
            case 6 % error correction
                corrError = bioinfoprivate.opttf(pval);
                if isempty(corrError)
                    error(message('bioinfo:seqlogo:InputOptionNotLogical', upper( char( okargs( k ) ) )));
                end
        end % end of switch
    end
end

%=============================================================
if isnumeric(P) %P has the profile
    if any(size(P,1)==[numAlphabet numAlphabet+1])
        isAA=true;
        freqM = P(1:numAlphabet, :)./PWMBg;
        freqM=freqM./sum(freqM);
        symbolList = cTypeChars(:);
    else
        isAA=false;
        if any(size(P,1)==[4 5])
            freqM=P(1:4, :);
        else
            error(message('bioinfo:seqlogo:IncorrectProfile'))
        end
    end
    % normalizing columns
    freqM = freqM./max(sum(freqM));
    corrError=false;
else % P has sequences that may include ambiguous symbols
    if iscell(P) || isfield(P,'Sequence')
        if isfield(P,'Sequence') % if struct put them in a cell
            P = {P(:).Sequence};
        end
        P = P(:);
        P = strrep(P,' ','-'); % padding spaces are not considered 'align' chars
        P = char(P); % now seqs must be a char array
    end

    if ~ischar(P)
        error(message('bioinfo:seqlogo:IncorrectInputType'))
    end
    seqs = upper(P);
    [numSeq, nPos] = size(seqs);
    
    % Create a profile count matrix  column: - positions
    % row:    - number of unique characters
    uniqueList = unique(seqs);
    
    if ~isAA && ~bioinfoprivate.isnt(uniqueList) && ~NTFlag
        warning(message('bioinfo:seqlogo:AmbiguousSequenceAlphabet'));
    end
    
    m = length(uniqueList);
    pcM = zeros(m, nPos);
    for i = 1:nPos
        for j = 1:m
            pcM(j,i) = sum(seqs(:,i) == uniqueList(j));
        end
    end
    
    % Compute the weight matrix used for graphically displaying the logo
    % Not considering wild card or gap YET, only for real symbols
    freqM = [];
    [symbolList, tmpIdx] = regexpi(uniqueList', '[A-Z]', 'match');
    symbolList = char(symbolList');

    if ~isempty(tmpIdx)
        for i = 1:length(tmpIdx)
            freqM(i, :) = pcM(tmpIdx(i),:); %#ok
        end
    end

    % The observed frequency of a symbol at a particular sequence position
    freqM = freqM/numSeq;
end

%===============================================================
%maxLen - the max sequence length in the set
maxLen = size(freqM, 2);

if endPos == 1
    endPos = maxLen;
end

% Check that startPos and endPos are not outside the seqs
if startPos > endPos
    error(message('bioinfo:seqlogo:StartAtGreaterThanEndAt', endPos));
end

if endPos > maxLen
    error(message('bioinfo:seqlogo:ExceedLimit', maxLen))
end

% == Compute the weight matrix used for graphically displaying the logo.
% Not considering wild card or gap YET, only for real symbols
freqM = freqM(:, startPos:endPos);
wtM = freqM; 
if isAA
    nSymbols = numAlphabet;
end

S_before = 0;
freqM(freqM == 0) = 1; % log2(1) = 0

% The uncertainty after the input at each position
S_after = -sum(log2(freqM./PWMBg).*freqM, 1);

if corrError
    % The number of sequences correction factor
    e_corr = (nSymbols -1)/(2* log(2) * numSeq);
    R = S_before - (S_after + e_corr);
else
    R = S_before - S_after;
end

R=R./max(R)*maxEntropy;

nPos = (endPos - startPos) + 1;
for i =1:nPos
    wtM(:, i) = wtM(:, i) * R(i);
end

if nargout >= 1
     % Create the seqLogo cell array
    W = cell(1,2);
    W(1,1) = {symbolList};
    W(1,2) = {wtM};
end

if isAA
    [wtM, symbolList] = sortWeightOrder(wtM, symbolList);
end

% Display logo
hfig = [];
if displayLogo
    wtM (wtM < 0) = 0;
    hfig = seqshowlogo(wtM, symbolList, isAA, startPos, cSeedStr);

    % if ~isempty(fileName)  % Save a image file
    %     hfig = seqshowlogo(wtM, symbolList, isAA, startPos, cSeedStr);
    % else
    %     hfig = seqshowlogo(wtM, symbolList, isAA, startPos);
    % 
    % end
end
end %seqlogo
%-------------------- Helper functions and callbacks ----------------%
function [p,s] = sortWeightOrder(weight, symbollist)
% Sort weight matrix by the sort of symbol list in ASCII direction order
% Here only needed for AA
[s, index] = sort(symbollist);
p=weight;
for i = 1:size(weight, 2)
    x=weight(:,i);
    p(:,i) = x(index);
end
end %sortWeightOrder
%--------------------------------------------------------------------%
function hFigure = seqshowlogo(varargin)
%SEQSHOWLOGO displays a Java seqlogo frame in a figure window
isAA = false;
seqType = 'NT';
filename = 'seqlogo.png'; %#ok!
saveLogo = false;%#ok!
wtMatrix = [];
symbols = [];
startPos = 1;

if nargin == 4 % Pass in weight Matrix, list of symbols and isAA
    wtMatrix = varargin{1};
    symbols = varargin{2};
    isAA = varargin{3};
    startPos = varargin{4}; 
elseif nargin == 5 % Pass in weight Matrix, list of symbols, isAA and filename
    saveLogo = true;%#ok!
    wtMatrix = varargin{1};
    symbols = varargin{2};
    isAA = varargin{3};
    startPos = varargin{4}; 
    cSeedStr = varargin{5};%#ok!


end

if isAA
    seqType = 'AA';
end

import com.mathworks.toolbox.bioinfo.sequence.*;
import com.mathworks.mwswing.MJScrollPane;
import java.awt.Dimension;
% Create the viewer
logoViewer = SequenceViewer(wtMatrix, symbols,startPos, seqType);
awtinvoke(logoViewer,'addSeqLogo()');
scrollpanel = MJScrollPane(logoViewer, MJScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,...
                              MJScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);

% Create a figure with the seqlogo panel on it and a uitoolbar
logoContainer = [];
hFigure = figure( ...
            'WindowStyle', 'normal',...
            'Resize', 'on', ...
            'Toolbar', 'none',...
            'NumberTitle','off',...
            'Tag', 'seqlogo',...
            'Name', cSeedStr,...
            'HandleVisibility', 'Callback',...
            'visible', 'off',...
            'DeleteFcn', {@onLogoClosing, logoViewer, logoContainer});

initFigureTools(hFigure, logoViewer)

% Set the figure widow size to fit the scrollPane
d = awtinvoke(scrollpanel, 'getPreferredSize()');
pos = getpixelposition(hFigure);
pos(3) = d.getWidth;
pos(4) = d.getHeight;
setpixelposition(hFigure,pos);
[logoP, logoContainer] = matlab.ui.internal.JavaMigrationTools.suppressedJavaComponent(scrollpanel, ...
                                       [0, 0, pos(3), pos(4)], hFigure);

set(logoContainer, 'units', 'normalized');
set(hFigure, 'visible', 'off')


end %seqshowlogo
%----------------------------------------------------------------------%
% % Using figure print function instead.
% % function printHandler(hsrc, event,logoViewer) %#ok
% % awtinvoke(logoViewer, 'logoPrint()');

%----------------------------------------------------------------------%
function saveHandler(hsrc, event, logoViewer) %#ok<INUSL>
awtinvoke(logoViewer, 'saveLogoDialog()')
end
%----------------------------------------------------------------------%
function onLogoClosing(hfig, event, logoViewer, logoContainer) %#ok<INUSL>
if ~isempty(logoViewer)
    awtinvoke(logoViewer, 'cleanup()');
    delete(logoContainer);
end
end
%--------------------------------------------------------------------
function initFigureTools(fig, logoViewer)
% helper function to set figure menus and toolbar
oldSH = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on')

% Handle toolbar 
toolbarHandle = uitoolbar('parent', fig);
hSave = uitoolfactory(toolbarHandle, 'Standard.SaveFigure');
set(hSave,  'ClickedCallback', {@saveHandler, logoViewer}, 'tooltip', 'Export Logo Image');

hPrint = uitoolfactory(toolbarHandle, 'Standard.PrintFigure');
set(hPrint, 'tooltip', 'Print');

% delete figure menus not used
%h1 = findall(fig,'Type','uimenu', 'Label','&Edit');
h1 = findall(fig,'Type','uimenu', 'Tag','figMenuEdit');
%h2 = findall(fig,'Type','uimenu', 'Label','&View');
h2 = findall(fig,'Type','uimenu', 'Tag','figMenuView');
%h3 = findall(fig,'Type','uimenu', 'Label','&Insert');
h3 = findall(fig,'Type','uimenu', 'Tag','figMenuInsert');
%h4 = findall(fig,'Type','uimenu', 'Label','&Tools');
h4 = findall(fig,'Type','uimenu','Tag','figMenuTools');
%h5 = findall(fig,'Type','uimenu', 'Label','&Desktop');
h5 = findall(fig,'Type','uimenu', 'Tag','figMenuDesktop');
delete([h1,h2,h3,h4,h5])

% Repair "File" menu
%hw = findall(fig,'Type','uimenu', 'Label','&File');
hw = findall(fig,'Type','uimenu', 'Tag','figMenuFile');
hf = get(hw,'children');
%h1 = findall(hw,'Label','&Save');
h1 = findall(hw,'Tag','figMenuFileSave');
%h2 = findall(hw,'Label','Print Pre&view...');
h2 = findall(hw,'Tag','figMenuFilePrintPreview');
%h3 = findall(hw,'Label','&Print...');
h3 = findall(hw,'Tag','printMenu');
%h4 = findall(hw,'Label', '&Close');
h4 = findall(hw,'Tag', 'figMenuFileClose');
delete(setxor(hf,[h1,h2,h3, h4]))

set(h1, 'label', '&Export Logo Image...', 'Callback', {@saveHandler, logoViewer});
set(h1,'Separator','on')

% Repair "Help" menu
%hw = findall(fig,'Type','uimenu','Label','&Help');
hw = findall(fig,'Type','uimenu','Tag','figMenuHelp');
delete(get(hw,'children'));
uimenu(hw,'Label','Bioinformatics Toolbox Help','Position',1,'Callback',...
       'helpview(fullfile(docroot,''toolbox'',''bioinfo'',''bioinfo.map''),''bioinfo_product_page'')')
uimenu(hw,'Label','Examples','Position',2,'Separator','on',...
       'Callback','demo(''toolbox'',''bioinfo'')')   
tlbx = ver('bioinfo');
mailstr = ['web(''mailto:bioinfofeedback@mathworks.com?subject=',...
           'Feedback%20for%20SeqLogo%20in%20Bioinformatics',...
           '%20Toolbox%20',tlbx(1).Version,''')'];
uimenu(hw,'Label','Send Feedback','Position',3,'Separator','on',...
       'Callback',mailstr);

set(0,'ShowHiddenHandles',oldSH)
end
