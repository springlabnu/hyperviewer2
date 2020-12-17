function [x, r, t] = NNLSCPU(A, b, CONFIG)
%NNLSCPU hyperspectral image cube decomposition via CPU implementation.
%   x = NNLSCPU(A, b) Returns a least-squares solution x to a set of linear 
%   equations Ax = b with the non-negativity constraint x >= 0. This 
%   function is a wrapper inteded for hyperspectral image decomposition
%   using NNLS. The FNNLS function is used if detected, which is about
%   2 times faster than the standard LSQNONNEG function. Otherwise
%   LSQNONNEG is used.
%  
%   x = NNLSCPU(A, b, CONFIG) uses the algorithm configuration specified by
%   CONFIG, 'fnnls' (ref 1) or the standard active-set method 'lsqnoneg' (ref 2)
%
%   [x, r, t] = NNLSCPU(...) also returns the residuals 
%   r = b - A*x where r is of size(b), and the total computation time t to
%   compute the solution x. 
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
%   [x, r, t] = NNLSCPU(A, b);
%
% Author: Eric Kercher, Northeastern University
% Date:   30JUN2020
%
% References
% [1] Bro, R., Jong, S. (1997). "A fast non‐negativity‐constrained least squares algorithm." J. Chemometrics.  11(5), 393-401.
% [2] Lawson and Hanson, "Solving Least Squares Problems", Prentice-Hall, 1974.

% Check arguments
if nargin < 2
    % Not enought args
    error('Not enough input arguments. A and b must be specified.');
    
elseif nargin > 3
    % Too many args
    error('Too many input agruments');
    
elseif nargin == 2
    % Set third arg to default value.
    try % Make sure fnnls exists on the matlab path
        fnnls(1,1);
        CONFIG = 'fnnls';
    catch % If not, default to lsqnonneg algorithm
        CONFIG = 'lsqnonneg';
    end
end

% System size
[n, m] = size(A);
if n < m
    error('System is underdetermined. n must be >= m.');
end

[imx, imy, num_wl] = size(b);
if n ~= num_wl
   error('Dimension mismatch. Number of spectral channels do not match.'); 
end
p = imx*imy;

% Reshape image cube into [n x p] matrix
h_b = reshape(b, [p n])';

% Check which pixels are zero valued
check_nonzero = sum(h_b, 1) > 0;

% Initialize
x_temp = zeros(m, p);
r_temp = zeros(n, p);

% Switch between two NNLS options
switch CONFIG
    case 'lsqnonneg'
        % Use LSQNONNEG for decomposition.
        tic
        parfor i = 1:p
            % Only run function if b is non-zero.
            if check_nonzero(i)
                % Run NNLS
                [x_temp(:, i), ~, r_temp(:, i), flag] = lsqnonneg(A, h_b(:,i));
                if ~flag
                    warning(['fit did not converge at pixel' i]);
                end
            end
        end
        t = toc;
        
    case 'fnnls'
        % Use FNNLS
        AtA = transpose(A) * A;
        Atb = transpose(A) * h_b;
        
        tic
        parfor i = 1:p
            % Only run function if b is non-zero.
            if check_nonzero(i)
                % Run FNNLS
                [x_temp(:, i), ~] = fnnls(AtA, Atb(:,i));
                
                % Collect result
                r_temp(:,i) = h_b(:,i) - A*x_temp(:, i);
            end
        end
        t = toc;
        
    otherwise
        
        error('Method not recognized. CONFIG must be ''lsqnonneg'' or ''fnnls''');
        
end

x = reshape(x_temp', [imx imy m]);
r = reshape(r_temp', [imx imy n]);

end