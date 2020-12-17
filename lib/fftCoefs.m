function T = fftCoefs(spec, n)

% Domain
lambda = linspace(0, 2*pi, size(spec, 1));

% initialize
T = zeros(length(n), size(spec, 2));

% Bandwidth
% w = 2*pi/(lambda(end) - lambda(1));

for i = 1:length(n)
    % Switch between sin and cos calculations
    if mod(i, 2) == 0
        T(i, :) = sin(n(i) * lambda) * spec ./ sum(spec, 1);
        
    else
        T(i, :) = cos(n(i) * lambda) * spec ./ sum(spec, 1);
        
    end
end
end