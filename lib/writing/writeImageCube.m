function [im, ME] = writeImageCube(im, filetype, savedir)

% Initailze error message.
ME = [];

cube = im.cube;
num_wl = size(cube, 3);

switch filetype
    case 'hyper' % Save hyperCFME Data
        
        a = 0; b = 0;
        % Loops backwards because hyperCFME image are saved in reverse
        % order
        for i = num_wl:-1:1

            filename = [im.label '-' num2str(im.power) 'mW_sys' ...
                         num2str(b) '_ch' num2str(a) '.tif'];
            
            image = uint16(cube(:, :, i) * 2^16);
            
            % Write
            imwrite(image, [savedir filename]);
            
            a = a + 1;
            if a == 4
                a = 0;
                b = b + 1;
            end
        end
        
        
end