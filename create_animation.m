%{
======================================================================
    Function that stitches together jpegs generated during the simulation into a .avi file.
    Adapted from https://www.mathworks.com/help/matlab/import_export/convert-between-image-sequences-and-video.html  
======================================================================
    INPUT:
        image_directory_location:    Location of directory containing images to be stitched together.
%}

function create_animation(image_directory_location)
    
    %% Find images in specified directory
    imageNames = dir(fullfile(image_directory_location, "*.jpg"));
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
    outputVideo = VideoWriter(fullfile(image_directory_location,'animation.mp4'), 'MPEG-4');
    
    % Specify frame rate
    outputVideo.FrameRate = 5; 
    
    % Open video, write frames, close video. 
    open(outputVideo)
    for i = 1:length(imageNames)
       img = imread(fullfile(image_directory_location,imageNames{i}));
       writeVideo(outputVideo,img)
    end
    close(outputVideo)

end