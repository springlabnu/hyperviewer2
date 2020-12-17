function [x3, r3, t] = phasorUnmix(S, p3, method)


% System size
[n, m] = size(S);
if n < m
    error('System is underdetermined. n must be >= m.');
end

[imx, imy, num_wl] = size(p3);
if n ~= num_wl
   error('Dimension mismatch. Number of spectral channels do not match.'); 
end

% Switch to 2D format
p2 = permute(p3,[3 2 1]);
p2 = p2(:,:);

% Calculate phasors

switch method
    case 'fourier'
        
        % Basis spectra define the anchor points in phasor space
        A = fftCoefs(S,  2:m);
        
        % Calculate each pixel phasor coordinates
        P = fftCoefs(p2, 2:m); 
        
        
    case 'chebyshev'
        
        % Basis spectra define the anchor points in phasor space
        A = fctCoefs(S,  2:m);
        
        % Calculate each pixel phasor coordinates
        P = fctCoefs(p2, 2:m); 
        
    otherwise
        error('Method not recognized');
end

x2 = zeros(m, imx*imy);

% Calculate umixing weights
tic
parfor i = 1:m
    
    % Logical indexing to loop through each spectra
    t = true(m, 1);
    t(i) = false;
    
    % Return weights as a simplex volume
    x2(i, :) = simplexvolume(A(:, t), P);
    
end
t = toc;

% Sum-to-one norm
x2 = x2 ./ repmat(sum(x2, 1), [m 1]);

% Residual
r2 = p2 - S*x2;


% Switch back to 3D image cube format
x3 = zeros(imx, imy, m);
r3 = zeros(imx, imy, n);
k = 1;
for i = 1:imx
    for j = 1:imy
        x3(i, j, :) = x2(:,k);
        r3(i, j, :) = r2(:,k);
        k = k + 1;
    end
end
end

