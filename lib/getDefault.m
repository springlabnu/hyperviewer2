function [out, flag] = getDefault(type, examp, varargin)

% Initialize
out = [];
flag = [];

switch type
    case 'demoDirectory'
        % Default directories for demo folder
        switch examp
            case 'hyper'
                out = [cd filesep 'demo' filesep 'mivic'];
                
            case 'dyes'
                out = [cd filesep 'demo' filesep 'soln'];
                
            case 'maestro'
                out = [cd filesep 'demo' filesep 'maestro'];
                
            case 'apples'
                out = [cd filesep 'demo' filesep 'apples'];
                
            case 'blood'
                out = [cd filesep 'demo' filesep 'BloodStain_RPI'];
                
            case 'oir'
                out = [cd filesep 'demo' filesep 'olympustest'];
                
            case 'oir-tiff'
                out = [cd filesep 'demo' filesep 'OlympusTiffs'];
                
            case 'czi'
                out = [cd filesep 'demo' filesep 'Zeiss'];
                
            otherwise
                flag = true;
                
        end
        
    case 'BSL'
        % Return path to default BSL mat file.
        out = fullfile('lib','BSL.mat');
        
    case 'spectra'
        % Default spectra to load for demos and filetypes
        % Return a cell array of strings listing the spectra names
        
        switch examp
            case '-'
                out = {};
                
            case 'hyper'
                % Default spectra to load if hyperCFME filetype detected
                out = {'AF633-7-14-20', 'AF660-7-14-20', ...
                    'AF680-7-14-20', 'AF700-v3-7-14-20', ...
                    'Laser10mW-7-14-20'};
                
            case 'dyes'
                out = {'AF633-20NOV', 'AF647-20NOV', ...
                    'AF660-20NOV', 'AF680-20NOV', ...
                    'AF700-20NOV', 'DarkPBS-110mW-20NOV'};
                
            case 'envi'
                
                out = [];
                
            case 'apples'
                
                out = {'Real Red', 'Real Green', 'Fake Red', 'Fake Green', ...
                    'WhiteRef_Apple'};
                
            case 'blood'
                % Default spectra to for bloodstain test case
                
                out = {'BloodStain', 'BlueStain', 'WhiteCloth', 'WhitePaper'};
                
            case 'maestro'
                % Default spectra to for maestro data
                out = {'BPD-TX', 'Bladder465', 'Organ465', 'Skin465',...
                    'Stomach465'};
                
            case 'oir'
                % Default spectra to for olympus data

                out = {'TSB-C1', 'TSB-C2', 'TSB-C3', 'TSB-C4', ...
                    'TSB-C5', 'OlympusBckgrd'};
                
            case 'oir-tiff'
                out = {'TSB-C1', 'TSB-C2', 'TSB-C3', 'TSB-C4', ...
                    'TSB-C5', 'OlympusBckgrd'};
                
            case 'macrotome_WL'
                out = {'brain_WL', 'OCT_WL', 'Sinus_WL', 'Skin_WL', ...
                    'Stomach_WL', 'Viable_WL', 'Necrosis_WL'};
                
            case 'macrotome_470'
                out = {'FITC_470', 'brain_470', 'OCT_470', 'skin_470', 'sinus_470', ...
                    'body_470', 'tumor_470', 'well_470'};
                
            case 'macrotome_635'
                % only 7 channels in cube so had to remove OCT_635 basis (max 7 bases)
                out = {'IR680_635', 'brain_635', 'Sinus_635', 'Skin_635', ...
                    'Stomach_635', 'Tumor_635', 'Well_635'};
                
            case 'czi-OPE'
                out = {'AF488-OPE-Zeiss', 'AF514-OPE-Zeiss', 'AF532-OPE-Zeiss', ...
                    'AF546-OPE-Zeiss', 'AF568-OPE-Zeiss', 'AF594-OPE-Zeiss', ...
                    'AF610-OPE-Zeiss', 'Background-OPE-Zeiss'};
                
            case 'czi-TPE'
                out = {'AF488-TPE-Zeiss', 'AF514-TPE-Zeiss', 'AF532-TPE-Zeiss', ...
                    'AF546-TPE-Zeiss', 'AF568-TPE-Zeiss', 'AF594-TPE-Zeiss', ...
                    'AF610-TPE-Zeiss', 'AF633-TPE-Zeiss', 'Background-TPE-Zeiss', ...
                    'BackgroundCell-TPE-Zeiss'};
                
            otherwise
                flag = true;
        end
        
        
    case 'roi'
        % Default ROI configurations
        % Return an roi object with mask
        
        imx = varargin{1}(1);
        imy = varargin{1}(2);
        
        switch examp
            % Change to load a default Polygon object
            case 'default'
                out.mask = ones(imx, imy);
                out.pts  = [0 0; 0 imx; imy imx; imy 0; 0 0];
                out.type = 'boundary';
                %out.axes = 'UICurrentAxes';
                
            case 'empty'
                out.mask = zeros(imx, imy);
                out.pts  = [];
                out.type = 'boundary';
                
            case 'fibercore'
                out.pts = [127.687500000000,115.046875000000;121.906250000000,168.234375000000;138.093750000000,224.890625000000;169.312500000000,259.578125000000;220.187500000000,286.171875000000;283.781250000000,287.328125000000;334.656250000000,267.671875000000;370.500000000000,227.203125000000;387.843750000000,168.234375000000;383.218750000000,124.296875000000;365.875000000000,84.9843750000000;318.468750000000,44.5156250000000;252.562500000000,27.1718750000000;195.906250000000,41.0468750000000;153.125000000000,72.2656250000000;127.687500000000,115.046875000000];
                out.mask = poly2mask(out.pts(:, 1), out.pts(:, 2), imx, imy);
                out.type = 'boundary';
                
            otherwise
                flag = true;
        end
    
    case 'roiColors'
        % Return a list of RGB triplets or color names in a cell array
        
        out = {[0.8 0.8 0.8], 'Green', 'Cyan', 'Yellow', 'Magenta', 'Red'};
        
    case 'samplingCorrection'
        % Return compression factor - scalar double. 
        switch examp
            case 'x'
                out = 1;
                
            case 'y'
                out = 6;
                
            otherwise 
                flag = true;
        end
        
    case 'phasorHarmonics'
        % Return phasor harmonic - scalar double
        switch examp
            case 'fourier'
                
                out = 1;
                
            case 'chebyshev'
                
                out = 2;
                
            otherwise
                flag = true;
        end
        
    case 'imageFilter'
        % Retun image filter computed by 'fspecial'
        switch examp
            
            case 'sigma'
                % Default width
                out = 4;
                
            case 'gaussian'
                % 6-by-6 square gaussian filter
                out = fspecial('gaussian', [6 6], varargin{1});
            % Add other image filters here
            otherwise
               flag = true;
        end
        
    case 'threshold'
        
        switch examp
            case 'default'
                out = 0;
            case 'hyper'
                out = 0;
            otherwise
                flag = true;
        end
        
    case 'parpool'
        
        switch examp
            case 'runOnStartup'
                out = 'Yes'; % Automatically start a parpool session
                % out = 'No'; % Do not start a parpool session
                % out = 'UserDecides'; % Ask user via popup box
        end
        
                
    otherwise
        % Type not found
        flag = true;
        
end
end