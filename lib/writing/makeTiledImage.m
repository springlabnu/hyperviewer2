function tiles = makeTiledImage(img, handles)
%% makeTiledImage(img, handles)

% Set default parameters
% TODO: Add inputdlg to let user define parameters
sizex = 182; % For hyperCFME - 182 pxls = 200 um
sizey = 182;
pxlbuffer = 0;
buffercolor = 'white';
direction = 'horizontal';

centerx = round((img.imx - sizex)/2);
centery = round((img.imy - sizey)/2);

% Default rectangle is centered
h = imrect(handles.axes1, [centery centerx sizex sizey]);
pos = wait(h);
delete(h);

x = [pos(1) pos(1)          (pos(1)+pos(3)) (pos(1)+pos(3)) pos(1)];
y = [pos(2) (pos(2)+pos(4)) (pos(2)+pos(4)) pos(2)          pos(2)];

mask = poly2mask(x, y, img.imx, img.imy);

% Find boundaries of mask
col1 = find(sum(mask, 1), 1, 'first');
col2 = find(sum(mask, 1), 1, 'last' );
row1 = find(sum(mask, 2), 1, 'first');
row2 = find(sum(mask, 2), 1, 'last' );

% Build tiles 

% Number of basis ims
m = handles.cfg.num_spec + 1;

% Calculate total size of tile
switch direction
    case 'horizontal'
        sizey       = sizey + pxlbuffer;
        sizex_total = sizex;
        sizey_total = sizey*m - pxlbuffer;       
    case 'vertical'
        sizex       = sizex + pixlebuffer;
        sizex_total = sizex*m - pxlbuffer;
        sizey_total = sizey;
        
end

%Initialize image as the buffer color
rgbvec = reshape(rgb(buffercolor), [1 1 3]);
tiles = repmat(rgbvec, [sizex_total, sizey_total]);

% Loop through images
for i = 1:m
    
    % Update new idx
    switch direction
        case 'vertical'
            x1 = (i - 1)*sizex + 1;
            x2 = x1 + sizex - pxlbuffer - 1;
            y1 = 1;
            y2 = sizey;
        case 'horizontal'
            x1 = 1;
            x2 = sizex;
            y1 = (i - 1)*sizey + 1;
            y2 = y1 + sizey - pxlbuffer - 1;
            
    end
    
    % Put composite first, then indv basis
    if i == 1
        % Compute composite image
        im = sum(cat(4, img.colorims.basis([handles.spec.filt]).im), 4);
        
        % Crop
        im = im(row1:row2, col1:col2, :);
        
        % Grab limts for contrast auto-adjust;
%         lims = stretchlim(im, 0);
    else
        im = img.colorims.basis(i - 1).im(row1:row2, col1:col2, :);
    end
    
    % Adjust contrast
%     im = imadjust(im, lims, []);
    
    % Place into tiled image
    tiles(x1:x2, y1:y2, :) = im;
    
    
end