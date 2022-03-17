% INSTRUCTIONS
% (1) Load data of interest into PeakSelector. Filter as appropriate.
%     (e.g., sigma rtNPh < 0.06, bounds on unwrapped z error, etc.)
% (2) Export data as an ASCII file with xy in pixels
% (3) Use custom tif to export the total raw data as a tiff file with
%     133.33 nm per pixel
% (4) Fill in appropriate filenames under USER PARAMETERS below for the
%     total raw data image and the ASCII file
%     (use the full path to the files)
% (5) Run script, which outputs a new ASCII file that can be re-loaded into
%     PeakSelector to create renderings without beads
%     (Import User ASCII, make sure to check box for headers, copy tab
%     deliminted 0:48 into the columns, xy coords in pixels)

clc, clear, close all % Start with a clean workspace
t.start = datetime('now'); % Measure how long pieces of this script take
% Can use between to measure times. To see total run time use:
% between(t.start,t.scriptFinished)

% Change Log
% RML 2022-02-08 clean up commenting

%% USER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add full paths to relevant files
totalRawFile = 'D:\rLee_localProjects\2022-03-15_Pekka_FilterFocalAdhesions\22.03.08-1\Run1-561\Run1-561_c123_sum_X14_processed_purged_IDL_totalRawData.tiff';
asciiFile = 'D:\rLee_localProjects\2022-03-15_Pekka_FilterFocalAdhesions\22.03.08-1\Run1-561\Run1-561_c123_sum_X14_processed_purged_IDL_ASCII.txt';

%%% Bead Finding
% These parameters are likely stable as long as the total raw tif stays
% 133.33 nm per pixel
rParticle = 4; % pixels, approximate size of beads
beadThresh = 10; % threshold on the bandpassed image

%%% Bead Removal
rRemoveX = 3; % pixels, how big of a region should be removed around each bead?
rRemoveY = 7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Loading

im = imread(totalRawFile,'tif');
im = max(im,[],3); % condense color image to grayscale
t.imageLoaded = datetime('now');

ascii = readmatrix(asciiFile);
headers = detectImportOptions(asciiFile);
headers = headers.VariableNames;
t.asciiLoaded = datetime('now');

%% Bead Finding
% bpass
im2 = padarray(double(im),10*rParticle*[1 1],0); % Padding allows finding beads at the edge
b = bpass(im2,0,rParticle); % input of 0 means no imaging smoothing
% imagesc(b)

% pknd (pixel-level)
pk = pkfnd(b,beadThresh,rParticle+1);

% cnt (sub-pixel)
cnt = cntrd(im2,pk,rParticle+11);
cnt(:,1) = cnt(:,1)-10*rParticle; % Remove the padding
cnt(:,2) = cnt(:,2)-10*rParticle;

%%%% Check peak finding
figure(1)
imshow(imadjust(im))
hold on
plot(cnt(:,1),cnt(:,2),'om','MarkerFaceColor','m','MarkerSize',4)
hold off
% title('Sanity Check: Are the beads identified?')

t.beadsFound = datetime('now');

%% Remove Bead Localizations

for kk = 1:length(cnt)

    x1 = cnt(kk,1);
    y1 = size(im,1)+1-cnt(kk,2); % Flip the coordinate system

    x = ascii(:,3);
    y = ascii(:,4);

    dx = x-x1;
    dy = y-y1;

    toRemove = ((dx.^2/rRemoveX.^2)+(dy.^2/rRemoveY.^2) < 1); % equation for an ellipse with width 2a & height 2b: x^2/a^2 + y^2/b^2 = 1

    ascii(toRemove,:) = [];

end

t.beadsRemoved = datetime('now');

%%%% Check Bead Removal
figure(2)
imshow(imadjust(im))
hold on
plot(ascii(:,3),size(im,1)+1-ascii(:,4),'.m','MarkerSize',2)
hold off
% title('Sanity Check: Were the beads removed?')

%% Save Work
T = array2table(ascii,'VariableNames',headers);
writetable(T,[asciiFile(1:end-4) '_beadsRemoved_rRemoveX' num2str(rRemoveX) '_rRemoveY' num2str(rRemoveY)  '.txt'],'delimiter','\t')

disp('Processing Completed')
t.scriptFinished = datetime('now');
