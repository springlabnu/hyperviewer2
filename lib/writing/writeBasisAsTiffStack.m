function writeBasisAsTiffStack(im, save_dir, subdirs, filename, app)
% WRITEBASISASTIFFSTACK   Tiff stack save function for HyperViewer2
%   writeBasisAsTiffStack(im, save_dir, subdirs, filename, app) Saves unmixed images to a Tiff stack
%       im              hyperim class object, unmixed by HyperViewer2
%       save_dir        directory to save image files
%       subdirs         1x3 cell of directory names for composite, individual color and individal grayscale image folders
%       filename        name of image to output - needs to be unique, overwriting is not checked
%       app             HyperViewer2 app instance

% Write colored composite basis-image

specfilt = 1:app.cfg.num_spec;
specfilt = specfilt(logical([app.spec.filt])); % In case filt was not a bool
colorcomp = im.colorImage(app, specfilt, 'x');
imwrite(colorcomp, [save_dir subdirs{1} filename 'basis_composite_image.tif']);

% Set up tiff objects for writing stacks

tColor = Tiff([save_dir subdirs{2} filename 'basis_stack_color.tif'],'w');
tBW = Tiff([save_dir subdirs{3} filename 'basis_stack_intensity_map.tif'],'w');

tagstructC.ImageLength = size(colorcomp, 1);
tagstructC.ImageWidth = size(colorcomp, 2);
tagstructC.Photometric = Tiff.Photometric.RGB;
tagstructC.BitsPerSample = 16;
tagstructC.SamplesPerPixel = 3;
tagstructC.Compression = Tiff.Compression.None;
tagstructC.SampleFormat = Tiff.SampleFormat.UInt;
tagstructC.RowsPerStrip = 512;
tagstructC.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstructC.ImageDescription = 'min=0.0 max=65535.0';

tagstructBW.ImageLength = size(colorcomp, 1);
tagstructBW.ImageWidth = size(colorcomp, 2);
tagstructBW.Photometric = Tiff.Photometric.MinIsBlack;
tagstructBW.BitsPerSample = 16;
tagstructBW.SamplesPerPixel = 1;
tagstructBW.Compression = Tiff.Compression.None;
tagstructBW.SampleFormat = Tiff.SampleFormat.UInt;
tagstructBW.RowsPerStrip = 512;
tagstructBW.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstructBW.ImageDescription = 'min=0.0 max=65535.0';

% Write tiff stacks, one slice per basis map

for i = 1:app.cfg.num_spec
    tColor.setTag(tagstructC);
    tBW.setTag(tagstructBW);
    write(tColor, uint16(im.colorImage(app, i, 'x') * (2^16-1)));
    write(tBW, uint16(im.x(:,:,i) * (2^16-1))), 
    if i ~= length(im.colorims.basis)
       tColor.writeDirectory();
       tBW.writeDirectory();
    end  
end

close(tColor);
close(tBW);