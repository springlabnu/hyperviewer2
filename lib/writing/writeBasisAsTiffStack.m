function writeBasisAsTiffStack(im, save_dir, subdirs, filename, handles)

% Write colored composite basis-image
specfilt = [handles.spec.filt];
try
imwrite( sum(cat(4, im.colorims.basis(specfilt).im), 4), ...
         [save_dir subdirs{1} filename 'basis_composite_image.tif'] );
catch
    imwrite( sum(cat(4, im.colorims.basis.im), 4), ...
         [save_dir subdirs{1} filename 'basis_composite_image.tif'] );
end

% Save spectral data in roi
% if true
%     save([save_dir 'analysis_data.mat'], 'im');
% end

tColor = Tiff([save_dir subdirs{2} filename 'basis_stack_color.tif'],'w');
tBW = Tiff([save_dir subdirs{3} filename 'basis_stack_intensity_map.tif'],'w');

tagstructC.ImageLength = size(im.colorims.basis(1).im, 1);
tagstructC.ImageWidth = size(im.colorims.basis(1).im, 2);
tagstructC.Photometric = Tiff.Photometric.RGB;
tagstructC.BitsPerSample = 16;
tagstructC.SamplesPerPixel = 3;
tagstructC.Compression = Tiff.Compression.None;
tagstructC.SampleFormat = Tiff.SampleFormat.UInt;
tagstructC.RowsPerStrip = 512;
tagstructC.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstructC.ImageDescription = 'min=0.0 max=65535.0';

tagstructBW.ImageLength = size(im.colorims.basis(1).im, 1);
tagstructBW.ImageWidth = size(im.colorims.basis(1).im, 2);
tagstructBW.Photometric = Tiff.Photometric.MinIsBlack;
tagstructBW.BitsPerSample = 16;
tagstructBW.SamplesPerPixel = 1;
tagstructBW.Compression = Tiff.Compression.None;
tagstructBW.SampleFormat = Tiff.SampleFormat.UInt;
tagstructBW.RowsPerStrip = 512;
tagstructBW.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstructBW.ImageDescription = 'min=0.0 max=65535.0';

for i = 1:handles.cfg.num_spec
    tColor.setTag(tagstructC);
    tBW.setTag(tagstructBW);
    write(tColor, uint16(im.colorims.basis(i).im * 65535));
    write(tBW, uint16(im.x(:,:,i) * 65535)), 
    if i ~= length(im.colorims.basis)
       tColor.writeDirectory();
       tBW.writeDirectory();
    end
     
end

close(tColor);
close(tBW);