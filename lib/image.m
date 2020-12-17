classdef image
    
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
        function flag = getImageCube(filetype)
            %% getImageCube -  Loads the hypersepctra image cube according to the
            % filetype and path provided in 'img'
            
            if image.isLoaded
                % Image already loaded
                return;
            end
            
            
            % Initailze error message.
            flag = [];
            
            % try % to load the image cube
            
            switch filetype
                case 'hyper' % Load hyperCFME Data
                    
                    % path specifies sub-directory containing images
                    if ~isfolder(image.path)
                        flag = 'HyperCFME data must be a directory';
                        return;
                    end
                    
                    % Get images in the sub directory
                    FileInfo = struct2cell(dir([image.path filesep '*mW_*.tif']));
                    if isempty(FileInfo)
                        flag = 'File not found';
                        image.filename = 'NOT_LOADED';
                        return;
                    end
                    
                    % Only need image names
                    image.imagefiles = FileInfo(1,:);
                    
                    % Number of lambda channels
                    num_wl = length(image.imagefiles);
                    
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
                    image.filename = [filenamebase ' - ' pwr_str ' mW'];
                    
                    % TODO: Test for non-integer powers
                    %     if length(pwr) > 4
                    %         if strcmp(pwr_str, 'point')
                    %             pwr(n) = str2double(pwr_str(6:end))*.1;
                    %         else
                    %             pwr(n) = str2double(pwr_str);
                    %         end
                    %     end
                    
                    % Load first image to get image size
                    I = double(imread([image.path filesep file1]));
                    [imx, imy] = size(I);
                    
                    % Initialize variables
                    im_cube   = zeros(imx, imy, num_wl);
                    L_sat_all = zeros(imx, imy);
                    
                    
                    for i = 1:num_wl
                        % Load Image
                        I = double(imread([image.path filesep imagefiles{i}]));
                        
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
                    filt = true(size(wlvec));
                    
                case 'maestro'
                    % path specifies sub-directory containing images
                    if ~isfolder(image.path)
                        error('Maestro data must be a directory');
                    end
                    
                    % Get images in the sub directory
                    FileInfo = struct2cell(dir([image.path filesep '*ms*.tif']));
                    
                    % Only need file names
                    imagefiles = FileInfo(1,:);
                    
                    num_wl = length(imagefiles);
                    
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
                    I = double(imread([image.path filesep file1]));
                    [imx, imy] = size(I);
                    
                    % Initialize variables
                    im_cube   = zeros(imx, imy, num_wl);
                    %         L_sat_all = zeros(imx, imy);
                    
                    wlvec = zeros(1, num_wl);
                    for i = 1:num_wl
                        % Path to image
                        impath = [image.path filesep imagefiles{i}];
                        
                        % Load Image
                        I = double(imread(impath));
                        
                        % Add the image to the hyperspectral image cube
                        im_cube(:,:,i) = I;
                        
                        % Parse filename
                        [~, filename, ~] = fileparts(impath);
                        parts = strsplit(filename, '_');
                        
                        % Save wavelength
                        wlvec(i) = str2double(parts{end});
                        
                        % TODO: Save exposure
                    end
                    
                    % Save base filename
                    filename = [sprintf('%s_', parts{1:end-2}), parts{end-1}];
                    
                    % TODO: Check for saturated pixels
                    SatPixels = ones(imx, imy);
                    
                    % No pwr to save
                    pwr = [];
                    filt = true(size(wlvec));
                    
                case 'envi' % Load ENVI data
                    % path specifies full path to the .dat file
                    % Requires use of enviread function
                    
                    [filepath, filename, ext] = fileparts(image.path);
                    if ~strcmp(ext, '.dat')
                        flag = 'Path to olympus datfile not specified';
                        return;
                    end
                    
                    hdrfile = [filepath filesep filename '.hdr'];
                    if exist(hdrfile, 'file') ~= 2
                        flag = 'Could not locate hdr file for specified dat file';
                        return;
                    end
                    
                    % Save image files
                    imagefiles = {[filename '.dat'], [filename '.hdr']};
                    
                    % Get image cube
                    [im_cube, info] = enviread(image.path, hdrfile);
                    
                    % For Specim images, rotate 90 degrees
                    if any(strcmp(info.interleave, {'BIL', 'bil'}))
                        im_cube = rot90(im_cube, -1);
                    end
                    
                    % Image parameters
                    [imx, imy, num_wl] = size(im_cube);
                    
                    % Parse sting of wavelengths into vector
                    wlvec = strsplit(info.wavelength, ',');
                    wlvec{1}(1) = []; wlvec{end}(end) = [];
                    wlvec = str2double(wlvec)';
                    if length(wlvec) ~= num_wl
                        % Some warning/error
                    end
                    
                    % For Trutag images, trim wl data from 460 - 700nm
                    if any(strcmp(info.interleave, 'bip'))
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
                    
                    
                    SatPixels = ones(imx, imy);
                    % No pwr to save
                    pwr = [];
                    filt = true(size(wlvec));
                    
                case 'oir'
                    % path specifies full path to oir file
                    % Requires use of the BioFormats Matlab Toolbox:
                    % https://www.openmicroscopy.org/bio-formats/downloads/
                    
                    [~, filename, ext] = fileparts(image.path);
                    if ~strcmp(ext, '.oir')
                        error('Path to olympus datafile not specified');
                    end
                    
                    imagefiles = {filename};
                    
                    % Load data
                    data = bfOpen3DVolume(image.path);
                    metadata = data{1, 2};
                    omeMeta  = data{1,4};
                    xmlString = char(omeMeta.dumpXML());
                    
                    % Save image cube
                    im_cube = data{1, 1}{1, 1};
                    [imx, imy, num_wl] = size(im_cube);
                    
                    % Get wavelength info
                    parts = strsplit(xmlString, '><');
                    
                    % Find correct tag
                    str = parts{startsWith(parts, 'ModuloAlongC')};
                    if isempty(str)
                        error('Could not find wl data in OME metadata');
                    end
                    
                    % Split string again to get info
                    info = strsplit(str, ' ');
                    
                    % Find wl Start
                    startstr = info{startsWith(info, 'Start="')};
                    temp = strsplit(startstr, '"');
                    wlstart = str2double(temp{2});
                    
                    % Find wl End
                    endstr = info{startsWith(info, 'End="')};
                    temp = strsplit(endstr, '"');
                    wlend = str2double(temp{2});
                    
                    % Find wl Step
                    stepstr = info{startsWith(info, 'Step="')};
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
                    filt = true(size(wlvec));
                    SatPixels = ones(imx, imy);
                    
                    % TODO: use this..
                    %         handles.oirSpectra = spectra;
                    
                case 'oir-tiff'
                    
                    % path specifies sub-directory containing images
                    if ~isfolder(image.path)
                        error('OIR-TIFF data must be a directory');
                    end
                    
                    % Get images in the sub directory
                    FileInfo = struct2cell(dir([image.path filesep '*lambda*.tif']));
                    
                    % Only need image names
                    imagefiles = FileInfo(1,:);
                    
                    % Number of lambda channels
                    % TODO: Check for extra files in folder
                    num_wl = length(imagefiles);
                    
                    % Take the first file to parse
                    file1 = FileInfo{1,1};
                    
                    % Get filename info, look for dash
                    parts = strsplit(file1, '-');
                    
                    % Save filename
                    filename = parts{1}(1:end-4);
                    
                    % Wavelength start/end surround the dash
                    wlstart = str2double(parts{1}(end-3:end));
                    wlend   = str2double(parts{2}(1:3)) - 20; % TODO: fix hardcode
                    
                    idx    = strfind(file1, 'nm step size');
                    wlstep = str2double(file1(idx-3:idx-1));
                    
                    wlvec = wlstart:wlstep:wlend;
                    
                    % Load first image to get image size
                    I = double(imread([path filesep file1]));
                    [imx, imy] = size(I);
                    
                    % Initialize variables
                    im_cube = zeros(imx, imy, num_wl);
                    
                    for i = 1:num_wl
                        % Path to image
                        impath = [path filesep imagefiles{i}];
                        
                        % Load Image
                        I = double(imread(impath));
                        
                        % Add the image to the hyperspectral image cube
                        im_cube(:,:,i) = I;
                    end
                    
                    % TODO: Check for saturated pixels
                    SatPixels = ones(imx, imy);
                    
                    % No pwr to save
                    pwr = [];
                    filt = true(size(wlvec));
                    
                case 'slice' % Load Davis Lab data
                    % path specifies full path to the .tif stack
                    
                    [filepath, filename, ext] = fileparts(image.path);
                    if ~strcmp(ext, '.tif')
                        flag = 'Path to .tif stack not specified';
                        return;
                    end
                    
                    % Save image files
                    imagefiles = {[filename '.tif']};
                    
                    % get wavelengths from .mat file in directory and put into vector
                    info = dir(filepath);
                    info(~endsWith({info.name},'.mat')) = [];
                    wlvec = struct2cell(load(fullfile(filepath, info(1).name)));
                    wlvec = wlvec{1};
                    wlvec = table2array(wlvec(:,1));
                    
                    % Image parameters
                    im_cube = imread(image.path,1);
                    [imx, imy] = size(im_cube);
                    num_wl = length(wlvec);
                    
                    % Get image cube
                    im_cube = ones(imx,imy,num_wl);
                    for i = 1:num_wl
                        im_cube(:,:,i) = imread(image.path,i);
                    end
                    
                    SatPixels = ones(imx, imy);
                    % No pwr to save
                    pwr = [];
                    filt = true(size(wlvec));
                    
                case 'czi'
                    % path specifies full path to czi file
                    % Requires use of the BioFormats Matlab Toolbox:
                    % https://www.openmicroscopy.org/bio-formats/downloads/
                    
                    [~, filename, ext] = fileparts(image.path);
                    if ~strcmp(ext, '.czi')
                        error('Path to zeiss datafile not specified');
                    end
                    
                    imagefiles = {filename};
                    
                    % Load data
                    data = bfOpen3DVolume(image.path);
                    metadata = data{1,2};
                    
                    % Save image cube
                    im_cube = data{1, 1}{1, 1};
                    [imx, imy, num_wl] = size(im_cube);
                    
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
                    filt = true(size(wlvec));
                    SatPixels = ones(imx, imy);
                    
            end
            
            % Initialize roi - default roi is entire image
            % image.roi = struct('mask',   ones(imx, imy), ... [dbl] binary mask, 1 inside the roi, 0 outside
            %                 'pts',    [0 0; 0 imx; imy imx; imy 0; 0 0], ... [dbl] (x,y) coords that define the roi boundary
            %                 'type',   'boundary');... [str] specify the type of roi {'boundary', 'point'}
            image.roi = getRegion([imx imy], 'default');
            
            image.roiStats = struct('savg',   [], ... [dbl] average spectrum in roi
                'sstd',   [], ... [dbl] standard dev of savg
                'xavg',   [], ... [dbl] average unmixing coefs in roi
                'xstd',   [], ... [dbl] standard dev of xavg
                'mchisq', [], ... [dbl] average reduced chi squared in roi
                'npx',    []);... [dbl] number of pixels in roi
                
            % OUTPUTS
            image.imx        = imx;
            image.imy        = imy;
            image.power      = pwr;
            image.wl         = wlvec;
            image.filename   = filename;
            image.imagefiles = imagefiles;
            image.cube       = im_cube;
            image.satpix     = SatPixels;
            image.filt       = filt;
            
            % Image is loaded
            image.isLoaded = true;
            image.isUnmixed = false;
            end
        end
        
    end