function T = phasorCoefs(spec, n, opt)


% Initialize
T = zeros(length(n), size(spec, 2));

switch opt
    
    case 'fourier'
        
        % Domain
        lambda = linspace(0, 2*pi, size(spec, 1));
        
        for i = 1:length(n)
            % Switch between sin and cos calculations
            if mod(i, 2) == 0
                T(i, :) = sin(n(i) * lambda) * spec ./ sum(spec, 1);
                
            else
                T(i, :) = cos(n(i) * lambda) * spec ./ sum(spec, 1);
                
            end
        end
        
    case 'chebyshev'
        % Domain
        lambda = linspace(-1, 1, size(spec, 1));
        
        for i = 1:length(n)
            
            T(i, :)  = (chebyshevU(n(i), lambda) .* sqrt(1 - lambda.^2)) * spec ./ sum(spec, 1);
            
        end
        
    otherwise
        error('Transform not recognized');

end
