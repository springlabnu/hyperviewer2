
if ~exist('im', 'var')
    [FileName1, PathName1] = uigetfile({'*.tif', 'Image Files (*.tif)'}, 'Choose the image to open.'); % Open Image File
    filename = [PathName1 FileName1];
    tstack = Tiff(filename,'r');

    [Width,Height] = size(tstack.read());
    k = length(imfinfo(filename));
    sz = imfinfo(filename).BitDepth;
    im = zeros(Width,Height,k);
    im(:,:,1)  = tstack.read();
    for n = 2:k
        tstack.nextDirectory()
        im(:,:,n) = tstack.read();
    end
    im = double(im/(2^sz-1)); % Convert to double image
end


results = unmix_accuracy(im);