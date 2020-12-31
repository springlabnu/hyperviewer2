function roi = getRegion(dims, opt, varargin)
%% getRegion - determine region mask and boundary points

% size of image
imx = dims(1);
imy = dims(2);

% Current axes
axes(gca);

switch opt
    
    case 'default'
        roi.mask = ones(imx, imy);
        roi.pts  = [0 0; 0 imx; imy imx; imy 0; 0 0];
        roi.type = 'boundary';
        
    case 'empty'
        roi.mask = zeros(imx, imy);
        roi.pts  = [];
        roi.type = 'boundary';
    
    case 'fibercore'
        roi.pts = [127.687500000000,115.046875000000;121.906250000000,168.234375000000;138.093750000000,224.890625000000;169.312500000000,259.578125000000;220.187500000000,286.171875000000;283.781250000000,287.328125000000;334.656250000000,267.671875000000;370.500000000000,227.203125000000;387.843750000000,168.234375000000;383.218750000000,124.296875000000;365.875000000000,84.9843750000000;318.468750000000,44.5156250000000;252.562500000000,27.1718750000000;195.906250000000,41.0468750000000;153.125000000000,72.2656250000000;127.687500000000,115.046875000000];
        roi.mask = poly2mask(roi.pts(:, 1), roi.pts(:, 2), imx, imy);
        roi.type = 'boundary';
        
    case 'point'
        % Get point
        [x, y] = ginput(1);
        roi.mask = zeros(imx, imy);
        roi.mask( round(y), round(x) ) = 1;
        roi.pts  = [round(x) round(y)];
        roi.type = 'boundary';
        
    case 'roi'
        % Get ROI
        [roi.mask, x, y] = roipoly;
        roi.pts = [x  y];
        roi.type = 'boundary';
        
    case 'scatter'
        % Get data from scatter plot
        pl = selectdata(varargin{1}{:});
        
        % Sometimes the ignore command doesn't work and 'pl' is returned as
        % a cell array. Check for empty entries and remove them
        % TODO: Hardcoded, fix this
        if iscell(pl)
%             idx = find(cellfun('isempty',pl));
%             pl(idx) = [];
            pl = pl{end};
        end
        
        % *** Note that mask and pts are created using the transpose
        % dimensions of the image. This is necessary because of the way
        % matlab indexes images vs. matricies (transposed).
        
        % Create pixel mask
        roi.mask = zeros(imy, imx);
        roi.mask(pl) = 1;
        roi.mask = roi.mask';
       
        % pixel points
        [x, y]   = ind2sub([imy imx], pl);
        roi.pts  = [x y];
        roi.type = 'points';
        
end