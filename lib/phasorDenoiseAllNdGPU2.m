function phasorDenoiseAllNdGPU2(app, type, n, sigma)
tic 

n_im = app.cfg.num_imgs;

[xdim,ydim,chan] = size(app.img(1).cube);
pixels = zeros(chan,xdim*ydim*n_im);
% number of pixels per image
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

% histogram phasor X,Y data
% values: array of bins, values of # of points that fall in that bin
% centers: x,y positions of centers of bins
[values, centers] = hist3([X' Y'], [numbins_x numbins_y]);

% create arrays of x positions, y positions
centers_X = centers{1,1};
centers_Y = centers{1,2};

% determine half the size of the bins
binsize_X = abs(centers_X(2) - centers_X(1)) / 2;
binsize_Y = abs(centers_Y(2) - centers_Y(1)) / 2;


% BIN PHASOR POINTS
[m,m2] = size(values);
% initialize array grouped_phasors to store the indexes of phasor points in
% a given bin set up by hist3 function
n_peakbin = max(max(values));
% grouped phasor indices (bin_x,bin_y,order added to list, phasor coords and phasor index)
grouped_phasors = zeros(m,m2,n_peakbin,n+1);
% also iteratively update the number of phasor pts in each bin to use as
% the index to add to grouped_phasors
num_phasors = zeros(m,m2);



for i = 1:size(X,2)
    % phasor coords ph_x and ph_y
    ph = T(:,i);
    ph_x = X(i);
    ph_y = Y(i);
    
    if isfinite(ph_x)
        % caclulate which bin the point is in
      
        bin_x = ceil((ph_x-centers_X(1)+binsize_X)/(2*binsize_X));
        bin_y = ceil((ph_y-centers_Y(1)+binsize_Y)/(2*binsize_Y));
       
        if bin_y > numbins_y
            bin_y = numbins_y;
        elseif bin_y == 0
            bin_y = 1;
        end
        
        if bin_x > numbins_x
            bin_x = numbins_x;
        elseif bin_x == 0
            bin_x = 1;
        end
        
        % add 1 to num_phasors in that bin and update list of phasor pt
        % indices in a given bin
        num_phasors(bin_x,bin_y) = num_phasors(bin_x,bin_y)+1;
        grouped_phasors(bin_x,bin_y,num_phasors(bin_x,bin_y),1:n) = ph;
        grouped_phasors(bin_x,bin_y,num_phasors(bin_x,bin_y),n+1) = i;
        
    end
end



% APPLY GAUSSIAN CONVOLUTION FILTER
% weights based on distance between pixels in phasor space
% sigma defines gaussian filter weights


for bin_x = 1:numbins_x
    
    %bin_x/numbins_x
    
    for bin_y = 1:numbins_y
        
        % create 3x3 array of bins to search (all points within 3*sigma of
        % current bin)
        
        bin_x_start = bin_x - 3;
        bin_x_end = bin_x + 3;
        cur_x = 4;
        if bin_x_start < 1
            bin_x_start = 1;
            cur_x = bin_x;
        end
        if bin_x_end > numbins_x
            bin_x_end = numbins_x;
        end
        
        bin_y_start = bin_y - 3;
        bin_y_end = bin_y + 3;
        cur_y = 4;
        if bin_y_start < 1
            bin_y_start = 1;
            cur_y = bin_y;
        end
        if bin_y_end > numbins_y
            bin_y_end = numbins_y;
        end
        
        num2search = num_phasors(bin_x_start:bin_x_end,bin_y_start:bin_y_end);
        bin_max = max(max(num2search));
        bins2search = grouped_phasors(bin_x_start:bin_x_end,bin_y_start:bin_y_end,1:bin_max,:);
        
     
        % get the values of the pixels corresponding to phasor pts in
        % bins2search
        s = size(num2search);
        px2search = zeros([s,bin_max,chan]);
        for search_x = 1:s(1)
            for search_y = 1:s(2)
                % consider each phasor point in current search bin
                for j = 1:num2search(search_x,search_y)
                    % index of phasor pt_i
                    pt_i = bins2search(search_x,search_y,j,n+1);
                    % pixel coords of pt
                    im_n = ceil(pt_i/n_pix);
                    px_x = ceil((pt_i-(im_n-1)*n_pix)/dim2);
                    px_y = mod(pt_i,dim2);
                    if mod(pt_i,dim2) == 0
                        px_y = dim2;
                    end

                    px2search(search_x,search_y,j,:) = app.img(im_n).disp(px_x,px_y,:);
                    
                end
            end
        end
                    
        bins2search = gpuArray(bins2search(:,:,:,1:n));
        
        % number of phasors in current bin to denoise
        n_bin = num_phasors(bin_x,bin_y);
        
        new_pixels = gpuArray(zeros(n_bin, chan));
        
        % new vectorized code
        px2search = gpuArray(px2search);
        for pt = 1:n_bin

            % phasor coords of pt stored in ph_pt
            ph_pt = squeeze(bins2search(cur_x,cur_y,pt,1:n));

            % for renormalization of current pixel
            norm_og = sum(px2search(cur_x,cur_y,pt,:));

            % phasor coords into matrix
            phasors = permute(bins2search,[4 3 2 1]);
            phasors = phasors(:,:);

            % pixel spectra values in same order as phasor coords
            pixels = permute(px2search,[4 3 2 1]);
            pixels = pixels(:,:);

            % calculate distance squared of this phasor from other phasors
            dist = sum((phasors(:,:)-ph_pt).^2);

            dist(dist > (9*sigma^2)) = NaN;

            weight = exp(-dist/(2*(sigma^2)));
            weight(isnan(weight)) = 0;
            norm = sum(pixels,1);
            norm(norm == 0) = 1;
            denoised_px = sum((pixels.*weight)./norm,2);

            new_pixels(pt,:) = denoised_px*norm_og/sum(denoised_px);
        end
        
        for pt = 1:n_bin
            % index of current phasor pt in the current bin
            pt_i = grouped_phasors(bin_x,bin_y,pt,n+1);

            % pixel coords of pt
            im_n = ceil(pt_i/n_pix);
            px_x = ceil((pt_i-(im_n-1)*n_pix)/dim2);
            px_y = mod(pt_i,dim2);
            if mod(pt_i,dim2) == 0
                px_y = dim2;
            end

            app.img(im_n).disp(px_x,px_y,:) = new_pixels(pt,:);
        end
    end
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