function denoised_cube = phasorDenoise(cube, X, Y)

%% denoise display image cube by phasor space nearest neighbor averaging
% inputs: raw image cube cube; phasor components X and Y as dbl vectors;
% density intensity values for each x,y point in heatmap as a dbl vector
% outputs: denoised cube denoised_cube with the same dimensions as the
% input cube

numbins = 50;

%{
% replace NaN values in X,Y with 0
for i = 1:length(X)
    if isnan(X(i))
        X(i) = 0;
    end
    if isnan(Y(i))
        Y(i) = 0;
    end
end
%}

% histogram phasor X,Y data
% values: array of bins, values of # of points that fall in that bin
% centers: x,y positions of centers of bins
[values, centers] = hist3([X' Y'], [numbins numbins]);

% create arrays of x positions, y positions
centers_X = centers{1,1};
centers_Y = centers{1,2};
% determins half the size of the bins
binsize_X = abs(centers_X(2) - centers_X(1)) / 2;
binsize_Y = abs(centers_Y(2) - centers_Y(1)) / 2;

%find bin with max number of phasor points and its center position
max = 0;
max_x = 0;
max_y = 0;
[m,n] = size(values);
for x = 1:m
    for y = 1:n
        if values(x,y) > max
            max = values(x,y);
            max_x = x;
            max_y = y;
        end
    end
end


% store the index of all phasor points that fall within a circle of radius
% radius centered at the center of the max bin

% circle radius is max of x,y binsizes
% might want to change this later
radius = min(binsize_X,binsize_Y);



% initialize arrays to store pixel locations corresponding to phasor points 
% within the histogram peak circle
peak_x = zeros(size(X));
peak_y = zeros(size(X));

% second dimension of input cube to be used to calculate pixel location
dim2 = size(cube,2);

%number of points in peak
n = 1;

for i = 1:size(X)
    
    % phasor coords ph_x and ph_y
    
    ph_x = X(i);
    ph_y = Y(i);
    
    % check that phasor pt is a number (not NaN)
    % check if distance to this point from center of max bin is less than
    % radius
    if (isfinite(ph_x) == 1) && ((ph_x-max_x)^2+(ph_y-max_y)^2)^(0.5) < radius)
        % convert from phasor index to pixel location
        x = idivide(i, int16(dim2)) + 1;
        y = mod(i,dim2);
        peak_x(n) = x;
        peak_y(n) = y;
        n = n+1;
    end
end

% now need to average 15 of these pixels
% sum each of the 16 channels of 15 of the pixels and average the values
peak_avg = cube(peak_x(1),peak_x(2),:);
for i = 2:15
    x = peak_x(i);
    y = peak_y(i);
    peak_avg = peak_avg + cube(x,y,:);
end
peak_avg = peak_avg/15;


% initilize a new denoised cube, same size as input cube
denoised_cube = zeros(size(cube));
% now reset all pixels in the phasor space peak circle to this peak_avg
% pixel value
for i = 1:n
    x = peak_x(i);
    y = peak_y(i);
    denoised_cube(x,y,:) = peak_avg;
end

% now average the rest of the pixels with their nearest neighbors in phasor
% space

for i = 1:size(X)
    % phasor coords ph_x and ph_y
    ph_x = X(i);
    ph_y = Y(i);
    
    % check that phasor pt is a number (not NaN)
    % check that this point isn't in peak bin
    if (isfinite(ph_x) == 1) && (((ph_x-max_x)^2+(ph_y-max_y)^2)^(0.5) >= radius)
        
        % convert from phasor index to pixel location
        x = idivide(i, int16(dim2)) + 1;
        y = mod(i,dim2);
        
        % sum all nearest neighbors pixel values and count number of
        % neighbors
        neighbors_sum = cube(x,y,:);
        num_neighbors = 1;
        
        % nearest neighbors are those within one radius of the point in
        % phasor space -- find nearest neighbors
        for j = 1:size(X)
            if ((ph_x-X(j))^2+(ph_y-Y(j))^2)^(0.5) < radius
                n_x = idivide(j, int16(dim2)) + 1;
                n_y = mod(j,dim2);
                neighbors_sum = neighbors_sum + cube(n_x,n_y,:);
                num_neighbors = num_neighbors + 1;
            end
        end
        
        % set pixel value in denoised cube to avg of nearest neighbors
        % pixels
        denoised_cube(x,y,:) = neighbors_sum/num_neighbors;
        
    end
end

    
        









%{
threshold_map = heatmap;

for i = 1:length(heatmap)
    if isnan(heatmap(i))
        heatmap(i) = 0;
    end
end


threshold = mean(heatmap);

for i = 1:length(heatmap)
    if heatmap(i)>threshold
        threshold_map(i) = 1;
    else
        threshold_map(i) = 0;
    end
end
%}






