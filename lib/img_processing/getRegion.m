function roi = getRegion(dims, opt, varargin)
%% getRegion - determine region mask and boundary points

% size of image
imx = dims(1);
imy = dims(2);

% Current axes
% axes(gca);
% TODO: need and option to specify which UIAxes was selected

switch opt
        
    case 'point'
        % Get point
        %[x, y] = ginput(1);
        %roi.mask = zeros(imx, imy);
        %roi.mask( round(y), round(x) ) = 1;
        %roi.pts  = [round(x) round(y)];
        axes = varargin{1};
        point = drawpoint(axes);
        roi.pts = point.Position;
        roi.mask = zeros(imx, imy);
        roi.mask( round(roi.pts(1)), round(roi.pts(2))) = 1;
        roi.type = 'boundary';
        
    case 'roi'
        % Get ROI
        % TODO: roipoly apparently doesn't work on UIAxes in AppDesigner,
        % update code. 
        % Update to drawpolygon
        %[roi.mask, x, y] = roipoly;
        
        %following hardcode of axis 1 needs to be changed at some point.
        %find a way to determine which image to put the ROI onto
        axes = varargin{1};
        poly = drawpolygon(axes);
        roi.pts = poly.Position;
        roi.mask = double(poly2mask(roi.pts(:,1),roi.pts(:,2),imx,imy));
        
        %Addition of initial point in order to have closed polygon plotted
        roi.pts(end+1,:) = roi.pts(1,:);
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
        
    otherwise
        
        % Try to load a default roi
        [roi, flag] = getDefault('roi', opt, dims);
        if flag
            warning('Default ROI not loaded properly.');
        end
end