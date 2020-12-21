function writeImages(im, save_dir, handles)

% Write colored composite raw-image
wlfilt = [im.filt];
imwrite( sum(cat(4, im.colorims.cube(wlfilt).im), 4), ... / sum(wlfilt), ...
         [save_dir 'raw_composite_image.tif'] );

% Write colored composite basis-image
specfilt = [handles.spec.filt];
try
imwrite( sum(cat(4, im.colorims.basis(specfilt).im), 4), ...
         [save_dir 'basis_composite_image.tif'] );
catch
    imwrite( sum(cat(4, im.colorims.basis.im), 4), ...
         [save_dir 'basis_composite_image.tif'] );
end
% % Write the saturated pixel mask image
% imwrite( uint16(L_sat_all .* 2^16), ...
%     [save_dir 'saturated pixel mask' '.tif'] );

% Write chi-squared map
f = figure('visible', 'off');
imagesc(im.redchisq);
axis square
set(gca, 'XTickLabel', [], 'YTickLabel', []);
set(f, 'Color', 'none');
colorbar; colormap jet;
export_fig(f, [save_dir 'ChiSquaredMap.png'], '-native');
clear f;

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

% Save spectral data in roi
if true
    save([save_dir 'analysis_data.mat'], 'im');
end

for i = 1:handles.cfg.num_spec
    
    % Save Color Images
    imwrite( im.colorims.basis(i).im * handles.state.brightness, ...
        [save_dir 'basis_image_' ...
        handles.spec(i).label '.tif'] );
    
    % Save Raw basis maps as 16-bit tiffs
    imwrite( uint16(im.x(:,:,i) * 2^16), ...
        [save_dir 'basis_image_intensity_map_' ...
        handles.spec(i).name '.tif'] );

end

if false
    % Generate tiled image
    tiles = makeTiledImage(im, handles);
    
    parts = strsplit(handles.expt.dirname, filesep);
    
    % Save image
    imwrite(tiles, [save_dir 'tiled_basis_images_' im.label '_' parts{end} '.tif']);
    
%     extra_dir = '/Users/eric/Documents/hyperCFME/data/5ColorDyeSoln-28ImgAvg-20NOV/Tiles/';
%     imwrite(tiles, [extra_dir 'tiled_basis_images_' im.label '_' parts{end} '.tif']);
end