function phasorDenoiseAllNdGPU(app, type, n, sigma)
tic 

n_im = app.cfg.num_imgs;

[xdim,ydim,chan] = size(app.img(1).cube);
pixels = zeros(chan,xdim*ydim*n_im);
n_pix = xdim*ydim;

for i = 1:n_im
    if ~app.img(i).isLoaded
        [app.img(i), ME] = getImageCube(app.img(i), app.expt.filetype);
        
        % Update filename before checking the error
        app.DirectoryListBox.Items = {app.img.filename};
        app.DirectoryListBox.ItemsData = 1:length({app.img.filename});

        if ~isempty(ME)
            app.msgbox.Text = ME;
            warning(ME);
        end
    end
    
    cube = app.img(i).cube;
    
    %going to put denoised cubes into app.img.disp
    app.img(i).disp = cube;
    
    %phantom images
    %convert uint16 to double
    cube = double(cube)/(2^16-1);
    cutoff = 0.013;
    pixmax = 1; %uint16
    
    %cutoff low intensity pixels (noise, no signal)
    mask = (max(cube,[],3)>cutoff);
    %mask = uint16(mask); %for uint16 cubes only
    cube2 = cube.*mask;

    %don't use saturated pixels
    mask2 = (max(cube,[],3) ~= pixmax);
    %mask2 = uint16(mask2); %for uint16 cubes only
    cube2 = cube2.*mask2;
    
    % get phasor coordinates
    % Reshape image and mask into columns of pixels
    im_pixels = permute(cube2, [3 2 1]);
    im_pixels = im_pixels(:,:);
    pixels(:,(i-1)*n_pix+1:i*n_pix) = im_pixels;
    
end

% if no signal, don't denoise
if max(max(pixels))==0
    return
end

if strcmp(type,'Fourier')
    T = fftCoefs(pixels, 1:n);
    %T = fftCoefs(pixels, 2:(n+1));
    
elseif strcmp(type,'Chebyshev')
    T = fctCoefs(pixels, 2:(n+1));
end

X = gpuArray(T(1,:));
Y = gpuArray(T(2,:));

% if there's only one data point to denoise, do nothing
if max(X) == min(X) && max(Y) == min(Y)
    return
end

% second (y) dimension of input cube to be used to calculate pixel location
dim2 = ydim;

% new array of denoised pixels
dn_pixels = pixels;

%make pixels a gpuArray
pixels = gpuArray(pixels);

for i = 1:size(X,2)
    x = X(i);
    y = Y(i);

    if isfinite(x)

        %get image numner and pixel coords corresponding to pixel i in
        %phasor list

        px = pixels(:,i);
        norm_og = max(px);
        
        % array D of distances from current pixel's phasor pt
        D = ((X-x).^2 + (Y-y).^2).^(0.5);
        % array of normalization factors
        N = sum(pixels,1);
        % array of weights for each pixel
        W = exp(-D/(2*sigma^2));
        W(isnan(W)) = 0;
        % array of normalized pixels
        norm_pixels = pixels ./ N;
        % recalculate pixel i's value by weighting normalized pixels
        denoised_px = sum(norm_pixels.*W,2);
        % reassign into dn_pixels array and renormalize
        dn_pixels(:,i) = denoised_px*norm_og/sum(denoised_px);
    end
end

% assign all dn_pixels into disp cube of associated image
for i = 1:size(X,2)
    im_n = ceil(i/n_pix);
    px_x = ceil((i-(im_n-1)*n_pix)/dim2);
    px_y = mod(i,dim2);
    if mod(i,dim2) == 0
        px_y = dim2;
    end
    app.img(im_n).disp(px_x,px_y,:) = dn_pixels(:,i);
end

%now unmix all images, skip pre-processing
for i = 1:n_im
    
     app.msgbox.Text = 'Unmixing all images...  0%';
            
    % Switch to axes3 to display unmixing stats
    hold(app.UIAxes3, "on");
    title(app.UIAxes3, 'Unmixing Speed');
    xlabel(app.UIAxes3, 'Image #');
    ylabel(app.UIAxes3, 'Time (s)');
    xlim(app.UIAxes3, [0 app.cfg.num_imgs]);
    legend(app.UIAxes3, 'off');
    set(app.UIAxes3, 'YScale', 'log');


    t_load = zeros(1, app.cfg.num_imgs);
    lag    = zeros(1, app.cfg.num_imgs);
    
    app.img(i).disp = im2double(app.img(i).disp);
    app.img(i) = app.img(i).getAmatrix(app);
    [app.img(i), ME] = unmixImage(app.img(i), app.cfg);
    if ~isempty(ME)
        app.msgbox.Text = app.msg.nogpu;
    end
    app.img(i) = app.img(i).postProcessImage(app);
    app.img(i).tmask = true(size(app.img(i).x));
    
    % Get region stats
    app.img(i).roiStats = getRegionStats(app.img(i), app.img(app.cfg.roiIdx).roi.mask);

    lag(i) = t_load(i) - app.img(i).evalTime;
    plot(app.UIAxes3, i, t_load(i), 'b*', i, lag(i), 'ro');
    ylim(app.UIAxes3, [-inf inf]);

    app.msgbox.Text = sprintf('Unmixing all images... %1.0f %%', ...
        app.expt.loadvec(i));
    drawnow;
end

hold(app.UIAxes3, "off");
app.msgbox.Text = 'Ready.';

toc