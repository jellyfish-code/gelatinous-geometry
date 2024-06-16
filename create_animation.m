%{
======================================================================
    Stitches together jpegs generated during the simulation into a .avi file.

    Adapted from https://www.mathworks.com/help/matlab/import_export/convert-between-image-sequences-and-video.html 
======================================================================
%}

% Specify data directory. Animation will be saved in the same directory.
% dataDir = fullfile([pwd,'/example_usage_data/offset_graft_elast0_0.05_elast1_0.05_vis_1000_bulk_mod_0.1']); 
dataDir = fullfile([pwd,'/Data/SLM_off_0319213413']); 


%% Find images
imageNames = dir(fullfile(dataDir, "*.jpg"));
imageNames = {imageNames.name}';

%% Sort image list in ascending order

% Extract number associated to each image
imageNumbers = split(imageNames, "."); 
imageNumbers = str2double(imageNumbers(:, 1)); 

% Sort based on number
[~, idx] = sort(imageNumbers);
imageNames = imageNames(idx);

%% Create animation
% Create VideoWriter object
outputVideo = VideoWriter(fullfile(dataDir,'animation.mp4'), 'MPEG-4');

% Specify frame rate
outputVideo.FrameRate = 5; 

% Open video, write frames, close video. 
open(outputVideo)
for i = 1:length(imageNames)
   img = imread(fullfile(dataDir,imageNames{i}));
   writeVideo(outputVideo,img)
end
close(outputVideo)