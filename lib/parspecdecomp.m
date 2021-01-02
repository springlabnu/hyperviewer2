function [x, r, t] = parspecdecomp(A, b, CONFIG, gpuCount)
%PARSPECDECOMP Parallel spectral decomposition of hyperspectral image
%cubes.
%   x = PARSPECDECOMP(A, b) decomposes the hyperspectral image cube b 
%   into a linear mixture of specified endmembers comprising A. A is the
%   endmember matrix containing m basis spectra in each column with n rows. 
%   b is the hyperspectral image cube of size [imx imy n]. The solution
%   x is a set of endmember maps of size [imx imy m] corresponding to each 
%   column of A. For each pixel, the least-squares solution x is computed  
%   for set of  linear equations Ax(i, j) = b(i, j) with the non-negativity 
%   constraint x >= 0. This function is a wrapper inteded for hyperspectral  
%   image decomposition using NNLS using parallel processing to accelerate 
%   computation.
%  
%   x = PARSPECDECOMP(A, b, CONFIG) uses the algorithm configuration 
%   specified by CONFIG:
%         'fnnls'     the fast-NNLS algorithm (ref 1)
%         'lsqnonneg' the standard active-set method (ref 2)
%         'gpu-fnnls' FNNLS implemented on a compatible GPU (ref 3)
%
%   x = PARSPECDECOMP(A, b, 'gpu-fnnls', gpuCount) deploys the NNLS
%   algorithm to a compatible CUDA-enabled GPU. To check for compatible
%   devices, run 'gpuDevice' in the command window. 
%
%   [x, r, t] = PARSPECDECOMP(...) also returns the residuals  r = b - A*x 
%   where r is of size(b), and the total computation time t to compute the
%   solution x in milliseconds. 
% 
%   Use with Parallel Computing Toolbox:
%       - Both LSQNONNEG and FNNLS functions are parfor compatible. A
%       significant speedup in eval time will occur if a parallel pool is
%       used. To start a parallel pool, type 'parpool' in the command
%       window. MATLAB will automatically start a parpool session when
%       'parfor' is called unless this feature is disables in the Parallel
%       Preferences settings. 
%
%   A = rand(8, 5);     % random endmember matrix
%   b = rand(6, 6, 8);  % random hyperspectral image cube
%   [x, r, t] = PARSPECDECOMP(A, b);
%
%   % Compare two algorithms
%   [x1, ~, t1] = parspecdecomp(A, b, 'fnnls');
%   [x2, ~, t2] = parspecdecomp(A, b, 'lsqnonneg');
%   t = t1/t2;                               % fnnls is t times faster than lsqnonneg
%   diff = sqrt(sum(sum(sum((x2 - x1).^2)))) % SSE between fnnls and lsqnonneg
%
% Author: Eric Kercher, Northeastern University
% Date:   30JUN2020
%
% References
% [1] Bro, R., Jong, S. (1997). "A fast non‐negativity‐constrained least squares algorithm." J. Chemometrics.  11(5), 393-401.
% [2] Lawson and Hanson, "Solving Least Squares Problems", Prentice-Hall, 1974.
% [3] Kercher, E, et al. "Video-rate hyperspectral unmixing for multiplexed molecular microscopy and microendoscopy". Nature Scientifc Reports (2020).

% Check arguments
if nargin < 2
    % Not enought args
    error('Not enough input arguments. A and b must be specified.');
    
elseif nargin > 4
    % Too many args
    error('Too many input agruments');
    
elseif nargin == 2
    % Default to CPU implementation first 
    CONFIG = 'fnnls';
    
elseif nargin == 3
    % If 'gpu-fnnls' is specified but not gpuCount, use the first GPU.
    % GPU count must be a vector.
    if strcmp(CONFIG, 'gpu-fnnls')
        if gpuDeviceCount == 0
            error('No compatible GPU on this machine.')
        end
        gpuCount = [1 0];
    end
    
elseif nargin == 4
    % Check for compatible GPU
    if strcmp(CONFIG, 'gpu-fnnls') && gpuDeviceCount == 0
        error('No compatible GPU on this machine.')
    end
    
end

% System size
[n, m] = size(A);
if n < m
    error('System is underdetermined. n must be >= m.');
end

[imx, imy, num_wl] = size(b);
if n ~= num_wl
   error('Dimension mismatch. A and b spectral resolution do not match.'); 
end
p = imx*imy;

% Reshape image cube into [n x p] matrix
h_b = reshape(b, [p n])';

% Switch between NNLS options
switch CONFIG
    case 'lsqnonneg'
        % Use LSQNONNEG for decomposition.
        % Check which pixels are zero valued
        check_nonzero = sum(h_b, 1) > 0;
        
        % Initialize
        x_temp = zeros(m, p);
        
        tic
        parfor i = 1:p
            % Only run function if b is non-zero.
            if check_nonzero(i)
                % Run NNLS
                [x_temp(:, i), ~, ~, flag] = lsqnonneg(A, h_b(:,i));
                if ~flag
                    warning(['fit did not converge at pixel' i]);
                end
            end
        end
        t = toc * 1000; % convert to ms
        
    case 'fnnls'
        % Use FNNLS
        
        % Check which pixels are zero valued
        check_nonzero = sum(h_b, 1) > 0;
        
        % Initialize
        x_temp = zeros(m, p);
        
        AtA = transpose(A) * A;
        Atb = transpose(A) * h_b;
        
        tic
        parfor i = 1:p
            % Only run function if b is non-zero.
            if check_nonzero(i)
                % Run FNNLS
                [x_temp(:, i), ~] = fnnls(AtA, Atb(:,i));
            end
        end
        t = toc * 1000; % convert to ms
        
    case 'gpu-fnnls'
        
        % Run kernel on gpu(s) using nnls mex file
        [x_temp, t] = nnls(single(A), single(A'*A)', single(h_b), int32(gpuCount));
        x_temp = double(x_temp);
        
    otherwise
        
        error('Method not recognized. CONFIG must be ''lsqnonneg'', ''fnnls'', or ''gpu-fnnls''');
        
end

% Compute residuals
r_temp = h_b - A*x_temp;

% Convert back to image cube dimensions
x = reshape(x_temp', [imx imy m]);
r = reshape(r_temp', [imx imy n]);

end

% FNNLS code included here for simplicity (from ref 1)
function [x, w] = fnnls(XtX, Xty)
tol = 10*eps*norm(XtX,1)*max(size(XtX));
[~,n] = size(XtX);
P = zeros(1,n);
z = zeros(n,1);
Z = 1:n;
x = P';
ZZ=Z;
w = Xty-XtX*x;

% set up iteration criterion
iter = 0;
itmax = 30*n;

% outer loop to put variables into set to hold positive coefficients
while any(Z) && any(w(ZZ) > tol)
    [~,t] = max(w(ZZ));
    t = ZZ(t);
    P(1,t) = t;
    Z(t) = 0;
    PP = find(P);
    ZZ = find(Z);
    nzz = size(ZZ);
    z(PP')=(Xty(PP)'/XtX(PP,PP)');
    z(ZZ) = zeros(nzz(2),nzz(1))';
    z=z(:);
    % inner loop to remove elements from the positive set which no longer belong
    
    while any((z(PP) <= tol)) && iter < itmax
        
        iter = iter + 1;
        QQ = find((z <= tol) & P');
        alpha = min(x(QQ)./(x(QQ) - z(QQ)));
        x = x + alpha*(z - x);
        ij = find(abs(x) < tol & P' ~= 0);
        Z(ij)=ij';
        P(ij)=zeros(1,max(size(ij)));
        PP = find(P);
        ZZ = find(Z);
        nzz = size(ZZ);
        z(PP)=(Xty(PP)'/XtX(PP,PP)');
        z(ZZ) = zeros(nzz(2),nzz(1));
        z=z(:);
    end
    
    x = z;
    w = Xty-XtX*x;
end
end