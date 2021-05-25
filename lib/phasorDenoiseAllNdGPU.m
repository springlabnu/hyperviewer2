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

X = T(1,:);
Y = T(2,:);

% if there's only one data point to denoise, do nothing
if max(X) == min(X) && max(Y) == min(Y)
    return
end


numbins_x = floor((max(X)-min(X))/(sigma));
numbins_y = floor((max(Y)-min(Y))/(sigma));

% second (y) dimension of input cube to be used to calculate pixel location
dim2 = ydim;

% new array of denoised pixels
dn_pixels = pixels;
% array J of index of phasor pt in X,Y arrays (make GPU array?)
J = (1:size(X,2));

for i = 1:size(X,2)
    x = X(i);
    y = Y(i);

    if isfinite(x)

        im_n = ceil(i/n_pix);
        px_x = ceil((i-(im_n-1)*n_pix)/dim2);
        px_y = mod(i,dim2);
        if mod(i,dim2) == 0
            px_y = dim2;
        end

        px = cube(px_x,px_y,:);
        norm_og = max(px);
        weight = 1;
        denoised_px = px*weight/norm_og;
        
        % array D of distances from current pixel's phasor pt
        D = ((X-x).^2 + (Y-y).^2).^(0.5);
        
        




% if numbins_x or numbins_y is 0 or 1, just apply filter to all points
% without binning and return
if numbins_x == 0 || numbins_x == 1 || numbins_y == 0 || numbins_y ==1
    for i = 1:size(X,2)
        x = X(i);
        y = Y(i);

        if isfinite(x)
            
            im_n = ceil(i/n_pix);
            px_x = ceil((i-(im_n-1)*n_pix)/dim2);
            px_y = mod(i,dim2);
            if mod(i,dim2) == 0
                px_y = dim2;
            end

            px = cube(px_x,px_y,:);
            norm_og = max(px);
            weight = 1;
            denoised_px = px*weight/norm_og;

            for j = 1:size(X,2)
                n_x = X(j);
                n_y = Y(j);

                % distance from current pixel
                d = ((x-n_x)^2+(y-n_y)^2)^(0.5);

                if d < 3*sigma
                    n_im_n = ceil(j/n_pix);
                    n_px_x = ceil((j-(im_n-1)*n_pix)/dim2);
                    n_px_y = mod(j,dim2);
                    if mod(j,dim2) == 0
                        n_px_y = dim2;
                    end

                    n_px = app.img(n_im_n).cube(n_px_x,n_px_y,:);
                    norm = max(n_px);
                    weight = exp(-(d^2)/(2*(sigma^2)));
                    denoised_px = denoised_px + n_px*weight/norm;
                end
            end
            app.img(im_n).disp(px_x,px_y,:) = denoised_px;
        end
    end
    return
end