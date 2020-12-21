function c_cube = smoothSpectra(wlvec, cube, specrange)
%% smoothSpectra - use a smoothing spline to smooth spectrum of pixels

[imx, imy, ~] = size(cube);
c_cube = zeros(size(cube));
wlvec = wlvec(:);

loadvec = linspace(0, 100, imx + 1);
loadvec(1) = [];
% handles.msgbox.String = sprintf('Calculating Splines ...  0%%');
drawnow;

% Default 4th order peicewise spline poly
X = [wlvec wlvec.^2 wlvec.^3 wlvec.^4];

% Per pixel smoothing spline fit
% NOTE: Curve fitting toolbox required for fit
for i = 1:imx
    parfor j = 1:imy
        
        % Get spectrum
        spec = squeeze(cube(i, j, :));
        
        % Fit
%         [f,gof,out] = fit(wlvec, spec, 'smoothingspline');
%         f = fit(wlvec, spec, 'smoothingspline');
        options = fitoptions('Method','Smooth','SmoothingParam',0.007);
        [f,~,~] = fit(wlvec, spec, 'smoothingspline', options);
%         out.p
        % Extract fit
        coefs = [f.p.coefs; f.p.coefs(end, :)];
        Y = sum(X.*coefs, 2);
        % Normalize to scale of spectra
        Y = Y .* max(spec) ./ max(Y);
        
        % Save fit into new cube
        c_cube(i, j, :) = Y;
        
    end
    
%     handles.msgbox.String = sprintf('Calculating Splines ...  %1.0f %%', ...
%                                     loadvec(i));
%     drawnow;

end