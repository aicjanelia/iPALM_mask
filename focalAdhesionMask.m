% INSTRUCTIONS
% (1) Load data of interest into PeakSelector. Filter as appropriate.
%     (e.g., sigma rtNPh < 0.06, bounds on unwrapped z error, etc.)
% (2) Export data as an ASCII file with xy in pixels
% (3) Optionally run beadRemoval_v0_anisotropic to remove beads
% (4) Fill in appropriate filenames under USER PARAMETERS below for the
%     focal adhesions image to use as a mask and the ASCII file
%     (use the full path to the files)
% (5) Run script, which saves a new ASCII file that can be re-loaded into
%     PeakSelector to create renderings without beads
%     (Import User ASCII, make sure to check box for headers, copy tab
%     deliminted 0:48 into the columns, xy coords in pixels)

clc, clear, close all % Start with a clean workspace
t.start = datetime('now'); % Measure how long pieces of this script take
% Can use `between` to measure times. To see total run time use:
% between(t.start,t.scriptFinished)

% Change Log
% RML 2022-03-15 focalAdhesionMask inspired by beadRemoval_v0_anisotropic

%% USER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add full paths to relevant files
mainDir = 'D:\rLee_localProjects\2022-03-15_Pekka_FilterFocalAdhesions\22.03.08-1\';
adhesionsFile = [mainDir 'pax\pax_c3.dat'];
asciiFile = [mainDir 'Run1-561\Run1-561_c123_sum_X14_processed_purged_IDL_ASCII_beadsRemoved_rRemoveX3_rRemoveY7.txt'];

% Thresholding
threshStrength = .8; % What fraction of an Otsu threshold to use --> lower is a more generous thresholding that will include more edge cases
minSize = 15; % Minimum size of a feature that is relevant (pixels^2)
dilateSize = 1; % Increasing the size of dilateSize will add extra regions around the focal adhesions that are included

%% Data Loading

% Read in the paxillin/adhesions image from a raw format
datSize = readmatrix([adhesionsFile(1:end-4) '.txt'],'NumHeaderLines',2);
datSize = [datSize(2) datSize(3)]; % Height, Width
fileID = fopen(adhesionsFile);
imThresh = fread(fileID,datSize,'uint16');
fclose(fileID);
imThresh = imrotate(imThresh,90);
t.imageLoaded = datetime('now');

% Read in the PeakSelector Data
ascii = readmatrix(asciiFile);
headers = detectImportOptions(asciiFile);
headers = headers.VariableNames;
t.asciiLoaded = datetime('now');

% Convert adhesions image to uint16 for processing
imThresh = (imThresh-min(imThresh(:)))./(max(imThresh(:))-min(imThresh(:))); % Scale 0 to 1
imThresh = uint16(imThresh*(2^16-1)); % 16 bit

%% Threshold the adhesions image

T = graythresh(imThresh); % Otsu threshold
BW = imbinarize(imThresh,threshStrength*T); % Apply threshold with a strength factor

BW = bwareaopen(BW, minSize); % Remove small objects
BW = imdilate(BW,strel('disk',dilateSize));

figure(1)
set(gcf,'Position',[400 200 800 1000]) % Get a good rough starting point for viewing this figure

subplot(2,2,1)
imshow(imadjust(imThresh))
title('Adhesions Image')

subplot(2,2,2) 
imshow(BW)
title('Mask')

subplot(2,2,3)
imshowpair(imadjust(imThresh),BW)
title('Overlay')

%% Only Keep adhesion Localizations

asciiIndex = round(ascii);
ind = sub2ind(size(BW),datSize(2)-asciiIndex(:,4),asciiIndex(:,3));
bwInd = BW(ind);

figure(1)
subplot(2,2,4)
imshow(imadjust(imThresh))
hold on
plot(ascii(bwInd,3),datSize(2)-ascii(bwInd,4),'.m','MarkerSize',1)
plot(ascii(~bwInd,3),datSize(2)-ascii(~bwInd,4),'.g','MarkerSize',1)
hold off
title('Localizations')

asciiKeep = ascii(bwInd,:);

%% Save Work
T = array2table(asciiKeep,'VariableNames',headers);
writetable(T,[asciiFile(1:end-4) '_masked.txt'],'delimiter','\t')

disp('Processing Completed')
t.scriptFinished = datetime('now');
