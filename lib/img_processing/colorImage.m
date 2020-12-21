function [colorIms, ME] = colorImage(im, rgbvec, normopt)
%% colorImage - convert intensity image into (imx x imy x 3) rgb images
% im is an image with k channels, colors is a (k x 3) matrix of rgb
% triplets specifying the color of the k'th image. 
[imx, imy, k] = size(im);

if size(rgbvec, 1) ~= k
	ME = 'Color map does not match image size.';
    return
end

% Initialize
colorIms = struct('im',    [], ... [dbl] (imx x imy x 3) rgb image
                  'rgb',   []);... [dbl] rgb trplet of colorname

% Convert color names to rgb triplets
rgbvec = reshape(rgbvec, [1 k 3]);

% Option to normalize image
if normopt
    % Normalize from 0 to 1
    im = mat2gray(im);
    % normIm = im / max(max(max(im)));
end

% Make pseudo-color image 
for i = 1:k
    
    % Multiply the rgb vector map by the image to pseudocolor
    colorIms(i).im = repmat(rgbvec(1, i, :), [imx imy]) .* ...
                     repmat(im(:, :, i), [1 1 3]);
                 
    colorIms(i).rgb = reshape(rgbvec(1, i, :), [1 3]);
    
end
% Options for normalization
% colorImTotal = colorImTotal / k;
% colorImTotal = mat2gray(colorImTotal);
% colorImTotal = colorImTotal /  max(max(max(colorImTotal)));