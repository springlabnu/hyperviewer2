function writeVideoFiles(sequence, vid_path, fps, handles)
% Wrapper for VideoWriter function

% Initialize
VidObj = VideoWriter( vid_path, 'MPEG-4');
VidObj.FrameRate = fps;
open(VidObj)

% Loop through image sequence, load each image, write it to the video.
for i = 1:length(sequence)
    
    % Skip empty frames
    if isempty(sequence{i})
        continue;
    end
    
    % Get frame and normalize
    frame = mat2gray(sequence{i} * handles.state.brightness);
    
    % Write frame to video
    writeVideo(VidObj, frame);
    
    % Inform user
    str = sprintf('Writing %d of %d images to %s...', ...
                  i, handles.cfg.num_imgs, vid_path);
    handles.msgbox.String = str;
    drawnow;
end
close(VidObj);