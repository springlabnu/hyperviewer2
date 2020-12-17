function T = fctCoefs(spec, n)

% Domain
lambda = linspace(-1, 1, size(spec, 1));

% Initialize
T = zeros(length(n), size(spec, 2));

for i = 1:length(n)
    
    T(i, :)  = (chebyshevU(n(i), lambda) .* sqrt(1 - lambda.^2)) * spec ./ sum(spec, 1);
    
end
end
