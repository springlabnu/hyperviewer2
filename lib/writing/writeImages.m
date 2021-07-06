function writeImages(im, save_dir, app)
% WRITEIMAGES   Save function for HyperViewer2
%   writeImages(im, save_dir, app) Saves unmixed image data to a directory for HyperViewer2
%       im              hyperim class object, unmixed by HyperViewer2
%       save_dir        directory to save image files
%       app             HyperViewer2 app instance

% Clear any files already in save_dir
delete([save_dir '*']);

% Write colored composite raw-image
wlfilt = [im.filt];
rawcomp = im.colorImage(app, wlfilt, 'cube');
imwrite(rawcomp,[save_dir 'raw_composite_image.tif'] );

% Write colored composite basis-image
specfilt = 1:app.cfg.num_spec;
specfilt = specfilt(logical([app.spec.filt])); % In case filt was not a bool
colorcomp = im.colorImage(app, specfilt, 'x');
imwrite(colorcomp, [save_dir 'basis_composite_image.tif']);

% Single channel images
for i = 1:app.cfg.num_spec
    % Save Color Images
    imwrite(im.colorImage(app, i, 'x'), [save_dir 'basis_image_' ...
        app.spec(i).label '.tif'] );
    
    % Save Raw basis maps as 16-bit tiffs
    imwrite(uint16(im.x(:,:,i) * (2^16-1)), [save_dir 'basis_image_intensity_map_' ...
        app.spec(i).name '.tif'] );
end

%% Optional additional files

% Save spectral data in roi as .mat file
if true
    save([save_dir 'analysis_data.mat'], 'im');
end

% Write chi-squared map
if false
    f = figure('visible', 'off');
    imagesc(im.redchisq);
    axis square
    set(gca, 'XTickLabel', [], 'YTickLabel', []);
    set(f, 'Color', 'none');
    colorbar; colormap jet;
    export_fig(f, [save_dir 'ChiSquaredMap.png'], '-native');
    clear f;
end

% Calculate diff between CPU & GPU
% if handles.gpuUnmix
if false
    
    x_gpu = im.x;  %#ok<*UNRCH>
    [x_cpu, ~, ~] = unmixSpeedTest(im.A, im.cube, 1, 'fnnls', [0 0]);
    
    % RMS Error
    rms = sqrt(mean((x_cpu{1} - x_gpu).^2, 3));
%     
%     % Write & save rmse map
    save([save_dir filesep 'RMSE_CPU_GPU.mat'], 'rms');
    
    f = figure('visible', 'off');
    imagesc(rms);
    axis image
    set(gca, 'XTickLabel', [], 'YTickLabel', []);
    set(f, 'Color', 'none');
    colorbar; colormap jet;
    export_fig(f, [save_dir 'CPU_GPU_difference.png'], '-native');
    clear f;
    
% else 
%     
%     x_cpu = handles.Coefcube{n};
%     x_gpu = [];
%     rms   = [];
%     
end

% Generate tiled image
if false
    
    tiles = makeTiledImage(im, app);
    
    parts = strsplit(app.expt.dirname, filesep);
    
    % Save image
    imwrite(tiles, [save_dir 'tiled_basis_images_' im.label '_' parts{end} '.tif']);
    
%     extra_dir = '/Users/eric/Documents/hyperCFME/data/5ColorDyeSoln-28ImgAvg-20NOV/Tiles/';
%     imwrite(tiles, [extra_dir 'tiled_basis_images_' im.label '_' parts{end} '.tif']);
end