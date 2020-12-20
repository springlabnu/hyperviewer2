classdef hyperim
    
    properties (Access = public)
       % Initialize image structure
        info  % [struct] directory info returend by dir
        path  % [str] path to image or image directory
        label % [str] name of image or image directory (displayed above image listbox)
        isLoaded % [logical] indicates if image has been loaded
        isUnmixed % [logical] indicates if image has been analyzed
        power % [dbl] laser power (if applicable)
        filename    % [str] image filename, no ext. (set as the string in the image listbox)
        imagefiles % [cell] full image filename(s) with extension
        cube % [dbl] raw hyperspectral data cube
        disp % [dbl] hyperspectral data cube after pre processing
        imx  % [dbl] image size in x dir
        imy  % [dbl] image size in y dir
        satpix % [logical] saturated pixel mask
        A % [dbl] unmixing matrix
        x % [dbl] unmixed basis image cube
        x_corr        % [dbl] unmixed basis image cube with flatfield correction
        ff_correction % [dbl] flatfield correction matrix
        redchisq % [dbl] map of reduced chi squared error in unmixing
        evalTime % [dbl] unmixing computation time
        colorims % [struct] pseudocolored composite images for display
        roi % [struct] roi mask and info
        roiStats  % [struct] statistics of the image defined by the roi
        wl % [dbl] vector of wavelength that define cube (can be different from basis spectra)
        phasor % [struct] phasor analysis variables
        filt  % [logical] vector used to filter wavlength channels if needed 
    end
    
    methods (Access = public)    
        function [obj, flag] = getImageCube(obj, filetype)
            %% getImageCube -  Loads the hypersepctra image cube according to the
            % filetype and path provided in 'img'
            
            if obj.isLoaded
                % Image already loaded
                return;
            end
            
            % Initialize error message.
            flag = [];
            
            % try % to load the image cube
            
            switch filetype
                case 'hyper' % Load hyperCFME Data
                    
                    % path specifies sub-directory containing images
                    if ~isfolder(obj.path)
                        flag = 'HyperCFME data must be a directory';
                        return;
                    end
                    
                    % Get images in the sub directory
                    FileInfo = struct2cell(dir([obj.path filesep '*mW_*.tif']));
                    if isempty(FileInfo)
                        flag = 'File not found';
                        obj.filename = 'NOT_LOADED';
                        return;
                    end
                    
                    % Only need obj names
                    obj.imagefiles = FileInfo(1,:);
                    
                    % Number of lambda channels
                    num_wl = length(obj.imagefiles);
                    
                    % Take the first file to parse
                    file1 = FileInfo{1,1};
                    
                    % Find the laser power in the filename
                    f = strfind(file1,'-');
                    if isempty(f)
                        f = strfind(file1, ' ');
                    end
                    
                    g = strfind(file1,'mW');
                    if isempty(g)
                        g = strfind(file1, 'Mw');
                    end
                    if isempty(g)
                        g = strfind(file1, 'MW');
                    end
                    % Save image file name
                    filenamebase = file1(1:f(end)-1);
                    
                    % Save laser power from filename
                    pwr_str = file1(f(end)+1:g(end)-1);
                    
                    % Assume integer powers for now, figure out non-integer
                    pwr = str2double(pwr_str);
                    
                    % Filename for listbox
                    obj.filename = [filenamebase ' - ' pwr_str ' mW'];
                    
                    % TODO: Test for non-integer powers
                    %     if length(pwr) > 4
                    %         if strcmp(pwr_str, 'point')
                    %             pwr(n) = str2double(pwr_str(6:end))*.1;
                    %         else
                    %             pwr(n) = str2double(pwr_str);
                    %         end
                    %     end
                    
                    % Load first image to get image size
                    I = double(imread([obj.path filesep file1]));
                    [obj.imx, obj.imy] = size(I);
                    
                    % Initialize variables
                    im_cube   = zeros(obj.imx, obj.imy, num_wl);
                    L_sat_all = zeros(obj.imx, obj.imy);
                    
                    
                    for i = 1:num_wl
                        % Load Image
                        I = double(imread([obj.path filesep obj.imagefiles{i}]));
                        
                        % Normalize Image (10 bit scale)
                        I = I / 1023;
                        
                        % Saturated pixels are assigned 1 in the current image
                        L_sat = (I > 1);
                        
                        % Running total of saturated pixels across all the channels
                        L_sat_all = L_sat_all + L_sat;
                        
                        % Add the image to the hyperspectral image cube
                        im_cube(:,:,i) = I;
                        
                    end
                    
                    % FLIP the image cube sequence and wavelengths to run from blue to red
                    im_cube = flip(im_cube, 3);
                    
                    % Normalize cube by laser power
                    %         im_cube = im_cube / pwr; % TODO: breaks when pwr is NaN
                    
                    % Get rid of hot corner
                    im_cube(347:end, 880:end, :) = 0;
                    
                    % Dilate and invert the saturation map such that saturated pixels are
                    % assigned 0. Make the saturation map binary
                    L_sat_all = (L_sat_all > 0);
                    
                    % Dilate the saturated pixels 4 pixels * 1.1 micron/pixel = 4 microns
                    L_sat_all = bwmorph(L_sat_all,'dilate', 4);
                    
                    % Now L_sat is assigned zeros for any saturated pixels across the image cube
                    SatPixels = (L_sat_all == 0);
                    
                    % No wavelength data embedded in images, hardcode
                    wlvec = [641 652 663 674 685 696 707 718 729 740 751 762 773 784 795 806]';
                    obj.filt = true(size(wlvec));
                    
                case 'maestro'
                    % path specifies sub-directory containing images
                    if ~isfolder(obj.path)
                        error('Maestro data must be a directory');
                    end
                    
                    % Get images in the sub directory
                    FileInfo = struct2cell(dir([obj.path filesep '*ms*.tif']));
                    
                    % Only need file names
                    obj.imagefiles = FileInfo(1,:);
                    
                    num_wl = length(obj.imagefiles);
                    
                    %         if num_wl ~= length(wlvec)
                    %             image_msg = sprintf(['Number of images: %d. \n'...
                    %                 'Expected number : %d. \n'...
                    %                 'Path: %s \n'], ...
                    %                 length(wlvec), num_wl,...
                    %                 [handles.Paths{n}]);
                    %             fprintf(image_msg);
                    %             error('Code did not find correct number of images to read');
                    %         end
                    %
                    % Take the first file to parse
                    file1 = FileInfo{1,1};
                    
                    % Load first image to get image size
                    I = double(imread([obj.path filesep file1]));
                    [obj.imx, obj.imy] = size(I);
                    
                    % Initialize variables
                    im_cube   = zeros(obj.imx, obj.imy, num_wl);
                    %         L_sat_all = zeros(imx, imy);
                    
                    wlvec = zeros(1, num_wl);
                    for i = 1:num_wl
                        % Path to image
                        impath = [obj.path filesep obj.imagefiles{i}];
                        
                        % Load Image
                        I = double(imread(impath));
                        
                        % Add the image to the hyperspectral image cube
                        im_cube(:,:,i) = I;
                        
                        % Parse filename
                        [~, obj.filename, ~] = fileparts(impath);
                        parts = strsplit(obj.filename, '_');
                        
                        % Save wavelength
                        wlvec(i) = str2double(parts{end});
                        
                        % TODO: Save exposure
                    end
                    
                    % Save base filename
                    obj.filename = [sprintf('%s_', parts{1:end-2}), parts{end-1}];
                    
                    % TODO: Check for saturated pixels
                    SatPixels = ones(obj.imx, obj.imy);
                    
                    % No pwr to save
                    pwr = [];
                    obj.filt = true(size(wlvec));
                    
                case 'envi' % Load ENVI data
                    % path specifies full path to the .dat file
                    % Requires use of enviread function
                    
                    [filepath, obj.filename, ext] = fileparts(obj.path);
                    if ~strcmp(ext, '.dat')
                        flag = 'Path to olympus datfile not specified';
                        return;
                    end
                    
                    hdrfile = [filepath filesep obj.filename '.hdr'];
                    if exist(hdrfile, 'file') ~= 2
                        flag = 'Could not locate hdr file for specified dat file';
                        return;
                    end
                    
                    % Save image files
                    obj.imagefiles = {[obj.filename '.dat'], [obj.filename '.hdr']};
                    
                    % Get image cube
                    [im_cube, obj.info] = enviread(obj.path, hdrfile);
                    
                    % For Specim images, rotate 90 degrees
                    if any(strcmp(obj.info.interleave, {'BIL', 'bil'}))
                        im_cube = rot90(im_cube, -1);
                    end
                    
                    % Image parameters
                    [obj.imx, obj.imy, num_wl] = size(im_cube);
                    
                    % Parse sting of wavelengths into vector
                    wlvec = strsplit(obj.info.wavelength, ',');
                    wlvec{1}(1) = []; wlvec{end}(end) = [];
                    wlvec = str2double(wlvec)';
                    if length(wlvec) ~= num_wl
                        % Some warning/error
                    end
                    
                    % For Trutag images, trim wl data from 460 - 700nm
                    if any(strcmp(obj.info.interleave, 'bip'))
                        idx1 = find(wlvec == 460);
                        if isempty(idx1)
                            idx1 = 1;
                        end
                        idx2 = find(wlvec == 700);
                        if isempty(idx2)
                            idx2 = length(wlvec);
                        end
                        
                        wlvec   = wlvec(idx1:idx2);
                        im_cube = im_cube(:, :, idx1:idx2);
                        num_wl  = length(wlvec);
                        
                    end
                    
                    
                    SatPixels = ones(obj.imx, obj.imy);
                    % No pwr to save
                    pwr = [];
                    obj.filt = true(size(wlvec));
                    
                case 'oir'
                    % path specifies full path to oir file
                    % Requires use of the BioFormats Matlab Toolbox:
                    % https://www.openmicroscopy.org/bio-formats/downloads/
                    
                    [~, obj.filename, ext] = fileparts(obj.path);
                    if ~strcmp(ext, '.oir')
                        error('Path to olympus datafile not specified');
                    end
                    
                    obj.imagefiles = {obj.filename};
                    
                    % Load data
                    data = bfOpen3DVolume(obj.path);
                    metadata = data{1, 2};
                    omeMeta  = data{1,4};
                    xmlString = char(omeMeta.dumpXML());
                    
                    % Save image cube
                    im_cube = data{1, 1}{1, 1};
                    [obj.imx, obj.imy, num_wl] = size(im_cube);
                    
                    % Get wavelength info
                    parts = strsplit(xmlString, '><');
                    
                    % Find correct tag
                    str = parts{startsWith(parts, 'ModuloAlongC')};
                    if isempty(str)
                        error('Could not find wl data in OME metadata');
                    end
                    
                    % Split string again to get info
                    obj.info = strsplit(str, ' ');
                    
                    % Find wl Start
                    startstr = obj.info{startsWith(obj.info, 'Start="')};
                    temp = strsplit(startstr, '"');
                    wlstart = str2double(temp{2});
                    
                    % Find wl End
                    endstr = obj.info{startsWith(obj.info, 'End="')};
                    temp = strsplit(endstr, '"');
                    wlend = str2double(temp{2});
                    
                    % Find wl Step
                    stepstr = obj.info{startsWith(obj.info, 'Step="')};
                    temp = strsplit(stepstr, '"');
                    wlstep = str2double(temp{2});
                    
                    % Save Wavelength
                    wlvec = wlstart:wlstep:wlend;
                    if num_wl ~= length(wlvec)
                        error('OIR Incorrect values found in metatdata');
                    end
                    
                    % Get dye info
                    k = 1;
                    chk = true;
                    spectra = cell(0);
                    while chk == true
                        str = ['dyeData name #' num2str(k)];
                        
                        if ~isempty(metadata.get(str))
                            spectra{k} = metadata.get(str);
                        else
                            chk = false;
                        end
                        k = k + 1;
                    end
                    
                    % Code to get all metadata
                    %         metadataKeys = metadata.keySet().iterator();
                    %         for i=1:metadata.size()
                    %             key = metadataKeys.nextElement();
                    %             value = metadata.get(key);
                    %             fprintf('%s = %s\n', key, value)
                    %         end
                    
                    % No pwr to save
                    pwr = [];
                    obj.filt = true(size(wlvec));
                    SatPixels = ones(obj.imx, obj.imy);
                    
                    % TODO: use this..
                    %         handles.oirSpectra = spectra;
                    
                case 'oir-tiff'
                    
                    % path specifies sub-directory containing images
                    if ~isfolder(obj.path)
                        error('OIR-TIFF data must be a directory');
                    end
                    
                    % Get images in the sub directory
                    FileInfo = struct2cell(dir([obj.path filesep '*lambda*.tif']));
                    
                    % Only need image names
                    obj.imagefiles = FileInfo(1,:);
                    
                    % Number of lambda channels
                    % TODO: Check for extra files in folder
                    num_wl = length(obj.imagefiles);
                    
                    % Take the first file to parse
                    file1 = FileInfo{1,1};
                    
                    % Get filename info, look for dash
                    parts = strsplit(file1, '-');
                    
                    % Save filename
                    obj.filename = parts{1}(1:end-4);
                    
                    % Wavelength start/end surround the dash
                    wlstart = str2double(parts{1}(end-3:end));
                    wlend   = str2double(parts{2}(1:3)) - 20; % TODO: fix hardcode
                    
                    idx    = strfind(file1, 'nm step size');
                    wlstep = str2double(file1(idx-3:idx-1));
                    
                    wlvec = wlstart:wlstep:wlend;
                    
                    % Load first image to get image size
                    I = double(imread([obj.path filesep file1]));
                    [obj.imx, obj.imy] = size(I);
                    
                    % Initialize variables
                    im_cube = zeros(obj.imx, obj.imy, num_wl);
                    
                    for i = 1:num_wl
                        % Path to image
                        impath = [obj.path filesep obj.imagefiles{i}];
                        
                        % Load Image
                        I = double(imread(impath));
                        
                        % Add the image to the hyperspectral image cube
                        im_cube(:,:,i) = I;
                    end
                    
                    % TODO: Check for saturated pixels
                    SatPixels = ones(obj.imx, obj.imy);
                    
                    % No pwr to save
                    pwr = [];
                    obj.filt = true(size(wlvec));
                    
                case 'slice' % Load Davis Lab data
                    % path specifies full path to the .tif stack
                    
                    [filepath, obj.filename, ext] = fileparts(obj.path);
                    if ~strcmp(ext, '.tif')
                        flag = 'Path to .tif stack not specified';
                        return;
                    end
                    
                    % Save image files
                    obj.imagefiles = {[obj.filename '.tif']};
                    
                    % get wavelengths from .mat file in directory and put into vector
                    obj.info = dir(filepath);
                    obj.info(~endsWith({obj.info.name},'.mat')) = [];
                    wlvec = struct2cell(load(fullfile(filepath, obj.info(1).name)));
                    wlvec = wlvec{1};
                    
                    % Image parameters
                    im_cube = imread(obj.path,1);
                    [obj.imx, obj.imy] = size(im_cube);
                    num_wl = length(wlvec);
                    
                    % Get image cube
                    im_cube = ones(obj.imx,obj.imy,num_wl);
                    for i = 1:num_wl
                        im_cube(:,:,i) = imread(obj.path,i);
                    end
                    
                    SatPixels = ones(obj.imx, obj.imy);
                    % No pwr to save
                    pwr = [];
                    obj.filt = true(size(wlvec));
                    
                case 'czi'
                    % path specifies full path to czi file
                    % Requires use of the BioFormats Matlab Toolbox:
                    % https://www.openmicroscopy.org/bio-formats/downloads/
                    
                    [~, obj.filename, ext] = fileparts(obj.path);
                    if ~strcmp(ext, '.czi')
                        error('Path to zeiss datafile not specified');
                    end
                    
                    obj.imagefiles = {obj.filename};
                    
                    % Load data
                    data = bfOpen3DVolume(obj.path);
                    metadata = data{1,2};
                    
                    % Save image cube
                    im_cube = data{1, 1}{1, 1};
                    [obj.imx, obj.imy, num_wl] = size(im_cube);
                    
                    % Get wavelength info - direct from channel "Names"
                    allKeys = arrayfun(@char, metadata.keySet.toArray, 'UniformOutput', false);
                    allValues = cellfun(@(x) metadata.get(x), allKeys, 'UniformOutput', false);
                    wlchrcell = allValues(contains(allKeys, 'Global Information|Image|Channel|Name'));
                    wlvec = flip(str2num(cell2mat(wlchrcell)));
                    
                    % Save Wavelength
                    
                    if num_wl ~= length(wlvec)
                        error('CZI Incorrect values found in metatdata');
                    end
                    
                    pwr = [];
                    obj.filt = true(size(wlvec));
                    SatPixels = ones(obj.imx, obj.imy);
                    
            end
            
            % Initialize roi - default roi is entire image
            % obj.roi = struct('mask',   ones(imx, imy), ... [dbl] binary mask, 1 inside the roi, 0 outside
            %                 'pts',    [0 0; 0 imx; imy imx; imy 0; 0 0], ... [dbl] (x,y) coords that define the roi boundary
            %                 'type',   'boundary');... [str] specify the type of roi {'boundary', 'point'}
            
            %obj.roi = getRegion([obj.imx obj.imy], 'default');
            
            obj.roiStats = struct('savg',   [], ... [dbl] average spectrum in roi
                'sstd',   [], ... [dbl] standard dev of savg
                'xavg',   [], ... [dbl] average unmixing coefs in roi
                'xstd',   [], ... [dbl] standard dev of xavg
                'mchisq', [], ... [dbl] average reduced chi squared in roi
                'npx',    []);... [dbl] number of pixels in roi
                
%             % OUTPUTS
%             obj.imx        = imx;
%             obj.imy        = imy;
             obj.power      = pwr;
             obj.wl         = wlvec;
%             obj.filename   = filename;
%             obj.imagefiles = imagefiles;
             obj.cube       = im_cube;
             obj.satpix     = SatPixels;
%             obj.filt       = filt;
%             
            % Image is loaded
            obj.isLoaded = true;
            obj.isUnmixed = false;
        end
        
        function obj = getAmatrix(obj, spec, filetype)
            %% getAmatrix - parses the basis spectra to match the image spectral
            % resolution (if required) and retuns d A matrix
    
            num_wl   = length(obj.wl);
            num_spec = length(spec);

            switch filetype
                case {'hyper', 'maestro'}

                    % Basis spectra already at correct resolution
                    obj.A = [spec.data];

                case {'envi', 'oir', 'oir-tiff', 'czi', 'slice'}

                    % Interpolate basis spec data to match image wl res
                    obj.A = zeros(num_wl, num_spec);
                    for i = 1:num_spec

                        % Interp the basis spectra at the image wavelength resolution
                        newspec = interp1(spec(i).wl, spec(i).data, obj.wl, 'pchip', 'extrap');

                        % Check if interpolation worked
                        if any(isinf(newspec)) || any(isnan(newspec))
                            % If not, warn user, A matrix kept at zero
                            warning('Ignoring spectra #%d; incorrect interpolation', i);

                        else
                            % All good
                            obj.A(:, i) = newspec;
                        end
                    end
            end

            % Normalize
            obj.A = obj.A ./ repmat(max(obj.A), [size(obj.A, 1), 1]);
        end
        
        function c_cube = smoothSpectra(obj, specrange)
            %% smoothSpectra - use a smoothing spline to smooth spectrum of pixels
            
            c_cube = zeros(size(obj.cube));
            wlvec = obj.wl(:);
            imgcube = obj.cube;
            
            % Default 4th order peicewise spline poly
            X = [wlvec wlvec.^2 wlvec.^3 wlvec.^4];
            
            % Per pixel smoothing spline fit
            % NOTE: Curve fitting toolbox required for fit
            for i = 1:obj.imx
                parfor j = 1:obj.imy
                    
                    % Get spectrum
                    spec = squeeze(imgcube(i, j, :));
                    
                    % Fit
                    %         [f,gof,out] = fit(wlvec, spec, 'smoothingspline');
                    %         f = fit(wlvec, spec, 'smoothingspline');
                    options = fitoptions('Method','Smooth','SmoothingParam',0.007);
                    [f,~,~] = fit(wlvec, spec, 'smoothingspline', options);
                    %         out.p
                    % Extract fit
                    coefs = [f.p.coefs; f.p.coefs(end, :)];
                    Y = sum(X.*coefs, 2);
                    % Normalize to scale of spectra
                    Y = Y .* max(spec) ./ max(Y);
                    
                    % Save fit into new cube
                    c_cube(i, j, :) = Y;
                    
                end
            end
        end
    end
end