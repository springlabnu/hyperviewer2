function roiStats = getRegionStats(im, mask)
%% getRegionStats - calculate means within an roi
% Get data for analysis
num_spec = size(im.x, 3);
num_wl   = size(im.disp, 3);

if isempty(im.cube)
    % Image not loaded
    return;
end

if im.imx ~= size(mask, 1) || im.imy ~= size(mask, 2) 
    % Image and mask dimension are not equal
    return;
end

% Number of pixels in ROI
num_px = sum(sum(mask));

% 3D masks
image_mask = logical(repmat(mask, [1 1 num_wl]));
basis_mask = logical(repmat(mask, [1 1 num_spec]));

% Keep spectral pixels included in mask (rows of pixels)
maskIm = reshape(im.disp(image_mask), [num_px num_wl]);

% Keep basis coeficients in mask (rows of coefs)
maskX = reshape(im.x(basis_mask), [num_px num_spec]);

if num_px == 1
    % Roi contatins one point, don't average
    savg = maskIm'; sstd = [];
    xavg = maskX';  xstd = [];
    
    % Normalize spectra and fit coefs (sum-to-one norm)
    smax = max(savg); %xsum = sum(xavg);
    savg = savg/smax; xavg = xavg/smax;
    
    avgRCS = im.redchisq(logical(mask));
    
else
    % Average pixel spectra in ROI
    savg = mean(maskIm)'; sstd = std(maskIm)';
    
    % Average x values in ROI
    xavg = mean(maskX)'; xstd = std(maskX)';
    
    % Normalize spectra and fit coefs
    smax = max(savg); %xsum = sum(xavg);
    savg = savg/smax; xavg = xavg/smax;
    sstd = sstd/smax; xstd = xstd/smax;
    
    % Average redchisq pixels
    avgRCS = mean(mean(im.redchisq(logical(mask))));
end

% OUTPUTS: Update handles structure
roiStats.savg   = savg;
roiStats.sstd   = sstd;
roiStats.xavg   = xavg;
roiStats.xstd   = xstd;
roiStats.mchisq = avgRCS;
roiStats.npx    = num_px;
