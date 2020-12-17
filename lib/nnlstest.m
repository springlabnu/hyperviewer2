function pass = nnlstest(dims, varargin)
% nnlstest is a function to quickly run the nnls mex file and test it's
% speed and accuracy. It creates a random data set with the specified width
% and number of pixels (or uses defaults if not given) and runs on both GPU
% and CPU and compares the result. If there is 
%
% 
%   Inputs:  dims  = a 3-element vector defining the hyperspectral
%                    parameters of the system:

%                    dims(1) = number of spectral channels (n)
%                    dims(2) = number of endmembers (m, m <= n)
%                    dims(3) = number of spatial pixels (p)

%            OPT   = String that turns on discriptions of the test. Default
%                    is 'verbose', which turns descriptions on. 'silent' 
%                    prints minimal information. 
%
%   Outputs: pass  = logical that is true if algorithm passes all tests.
%                    Otherwise it is false.
%
% Author: Eric Kercher, Northeastern University
% Date:   30AUG2017

%% Parse inputs
if nargin == 0 
    % Run all defaults
    dims = [20 20 gpuDeviceCount];
    OPT1 = 'default';
    OPT2 = 'CPU_on';
    
elseif nargin == 1
    % Use default system sizes and assign OPT
%     dims = [20 20 1];
    OPT1 = 'default';
    OPT2 = 'CPU_on';
    
elseif nargin == 2
    % Parameters are specified, assign OPT as defaults
    OPT1 = varargin{1};
    OPT2 = 'CPU_on';
    
elseif nargin == 3
    
    OPT1 = varargin{1};
    OPT2 = varargin{2};
    
elseif nargin > 4
    % Too many args
    error('Too many input arguments.');
    
end

switch OPT1
    case 'default'
        % Set default value of verbose logical
        verbose = true;
    case 'verbose'
        verbose = true;
    case 'silent'
        verbose = false;
    otherwise 
        verbose = false;
end

switch OPT2
    case {'default', 'CPU_on'}
        cpucomp = true;

    case 'CPU_off'
        cpucomp = false;

    otherwise 
        cpucomp = true;
end


% Logical to keep track of tests
pass = true;

% Define parameters
% dims = squeeze(dims);
if length(dims) == 3
    n = dims(1);
    m = dims(2);
    p = dims(3);
else
    error(['Missing dimensions. You must specify 3 hyperspectral ' ...
        'dimensions to define the system.']);
end

if m > n
    error(['System is underdetermined; the number of spectral channels '...
          'must be greater than or equal to the number of endmembers']);
end



% First check for compatible GPU(s)
if gpuDeviceCount < 0    
    pass = false;
    error('No compatible GPUs were found on this machine.');
elseif (gpuDeviceCount > 0) && verbose
    fprintf('GPU''s found: %d \n', gpuDeviceCount);
end

if verbose
    fprintf('System size: %d channels, %d endmembers\n', n, m);
    fprintf('# of pixels: %d \n\n', p);
    fprintf('Generating random system...');
end

% Create random data from system parameters
A   = rand(n, m);
AtA = A'*A;
b   = rand(n, p);
if verbose
    fprintf(' Done. \n');
end


% Does it run?
idx = zeros(gpuDeviceCount, 1);
idx(1) = 1;
try
    [x_gpu, ~] = nnls(single(A), single(AtA)', single(b), int32([1 0]));
    if verbose
        fprintf('NNLS kernel ran successfully. ');
    end
catch ME
    warning('NNLS ran unsuccessfully, nnls exited with the following error:');
    throw(ME);
end


%% Test 1: Accuracy vs. CPU

if cpucomp

if verbose
    fprintf('Comparing results to CPU: \n\n');
end

% Run on CPU
% tic
% x_cpu = cell2mat(cellfun(@(y) lsqnonneg(A, y), num2cell(b, 1), 'UniformOutput', false))';
% t_cpu = toc;
% x_cpu = x_cpu';

[x_cpu, ~, t_cpu] = NNLSCPU(A, permute(b, [3 2 1]), 'fnnls');
x_cpu = permute(x_cpu, [3 2 1]);

% Define tolerance
tol = 1e-6;
if verbose
    fprintf('Tolerance:         %3.2e \n', tol)
end

% Calculate norm: an average of the difference between each pixel, divided
% by the number of pixels
norm = mean(sum(abs(x_gpu - x_cpu), 2)) / p;

if verbose
    fprintf('Avg Pixel Norm:    %3.2e \n', norm);
end

if norm > tol 
    warning('System check failed. GPU result is likely incorrect.');
    pass = false;
elseif (norm < tol) && verbose
    fprintf('System check passed. \n\n');
end

else
   if verbose
       fprintf('CPU Comparison turned off. \n\n');
   end
end

%% Test 2: Run time 
% if gpuDeviceCount > 1
    if verbose
        fprintf('Running GPU speed test. Decomp Stats: \n');
        fprintf(' # GPU''s     Time (ms)      Speed (fps)   Speedup Factor \n');
        fprintf(' --------   -----------    -----------   -------------- \n\n');
    end
    
    nrep = 10;
    t = zeros(nrep, 1);
    gpuFilt = int32(zeros(gpuDeviceCount, 1));
    for i = 1:gpuDeviceCount
        for j = 1:nrep
            gpuFilt = int32(zeros(2, 1));
            gpuFilt(i) = int32(1);
            
            [x_gpu_c, t(j)] = nnls(single(A), single(AtA)', single(b), gpuFilt);
            
            % Check that results agree with original run
            if ~isequal(x_gpu, x_gpu_c)
                pass = false;
            end
        end
        t_avg   = mean(t);
        t_std   = std(t);
        fps     = 1000/t_avg;
        
        if cpucomp
        xfactor = t_cpu*1000/t_avg;
        
        if verbose
            fprintf('     %d     %3.1f +/- %2.1f       %3.1f            %4.0f \n', ...
                i, t_avg, t_std, fps, xfactor);
            
        else
            fprintf(' %d:  %3.1f +/- %2.1f \n', i, t_avg, t_std);
        end
        
        else 
            if verbose
            fprintf('     %d     %3.1f +/- %2.1f       %3.1f            N/A \n', ...
                i, t_avg, t_std, fps);
            
        else
            fprintf(' %d:  %3.1f +/- %2.1f \n', i, t_avg, t_std);
            end
        end
            
    end
    
    
%     % Run through backwards
%     for i = gpuDeviceCount:-1:1
%         x_gpu_c = nnls(single(A)', single(AtA)', single(b), i);
%         % Check that results agree with original run
%         if ~isequal(x_gpu, x_gpu_c)
%             pass = false;
%         end
%     end
%     
% end

if pass && verbose
    fprintf(' \nAlgorithm passed all nnlsmexTest checks. \n');
elseif ~pass && verbose
   fprintf(' \nAlgorithm did not pass all nnlsmexTest checks. \n'); 
end
