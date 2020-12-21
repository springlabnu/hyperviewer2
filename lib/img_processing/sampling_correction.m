function c_cube = sampling_correction(cube, p)
%% sampling_correction - 
% Compresss an image by averaging px pixels in the x and py pixels in the
% y dimensions for z replicates

isLogic = islogical(cube);
cube    = double(cube);
[imx, imy, imz] = size(cube);

% Compressed image dimensions
imx = floor(imx / p.x);
imy = floor(imy / p.y);
c_cube = zeros(imx, imy, imz);


for i = 1:imz
    temp = conv2(cube(:,:,i), ones(p.x, p.y), 'valid');
    c_cube(:,:,i) = temp(1:p.x:end, 1:p.y:end) / p.x*p.y;
end

% Reset to logical
if isLogic
    c_cube = logical(c_cube);
end