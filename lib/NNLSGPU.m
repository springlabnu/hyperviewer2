function [x, r, t] = NNLSGPU(A, b, gpuCount)
%NNLSGPU Returns a least-squares solution to a set of linear of equations 
% Ax = b with the non-negativity constraint x >= 0. This function is
% written for hyperspectral image decomposition. NNLSGPU calls nnls.mex64 
% file which loads and runs the nnls kernel on a CUDA-enabled GPU.
%
%  INPUTS:  A        = n by m matrix
%           b        = hyperspectral image cube (imx x imy x n)
%           gpuCount = (Optional) Integer specifying the number of GPUs to
%                      use for the NNLS algorithm. Default is the total
%                      GPU count on the system.
% 
%  OUTPUTS: x = hyperspectral solution cube (imx x imy x m)
%           r = residual norm of solution, norm = Ax - b (imx x imy x m)
%           t = Kernel evlauation time in seconds
%
% Author: Eric Kercher, Northeastern University
% Date:   30JUN2020
%

% Check GPU status
if gpuDeviceCount == 0
    error('Could not find compatible GPU on this machine')
end

% Check arguments
if nargin < 2
    % Not enought args
    error('Not enough input arguments. A and b must be specified.');
    
elseif nargin > 3
    % Too many args
    error('Too many input agruments');
    
elseif nargin == 2
    % Set third arg to default value.
    gpuCount = gpuDeviceCount;
end

% Get system parameters
[n, m] = size(A);
[imx, imy, ~] = size(b);

% Reshape image cube into [n x p] matrix
b_in = reshape(b, [imx*imy n])';

% Run kernel on gpu(s) using mex file
[x_gpu, t] = nnls(single(A), single(A'*A)', single(b_in), int32(gpuCount));

% Compile results
x = reshape(double(x_gpu)', [imx imy m]);
r = reshape(double(b_in - A*x_gpu)', [imx imy n]);
end
