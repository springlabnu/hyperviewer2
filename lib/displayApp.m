function displayApp(app, varargin)
%% displayApp(app) - update the entire hyperviewer2 GUI display
% displayApp(app, 'ax1', ..., 'ax3') allows user to select which UIAxes
% to update. Some actions only affect 1 or 2 UIAxes, this allows the gui to
% skip unchanged dispay objects to increase speed. For example: 
% 
% displayApp(app, 'ax1') only updates UIAxes1
% 
% displayApp(app, 'ax2', 'ax3') updates UIAxes2 and UIAxes2
%
% displayApp(app, 'ax1, 'ax2', 'ax3') is the same as displayApp(app)

% If no experiment, do nothing
if ~app.expt.loaded
    app.msgbox.Text = app.msg.loadexpt;
    return;
end

% Current selected image and spectra
n = app.state.img;
m = app.state.spec;

if strcmp(app.img(n).filename, 'NOT_LOADED')
    % File was not loaded properly the first time, skip it
    app.msgbox.Text = 'File not loaded properly.';
    return;
end

% Check if image has been loaded
if ~app.img(n).isLoaded
    % If not, load image
    
    % Inform user
    app.msgbox.Text = app.msg.loadingcube;
    
    [app.img(n), ME] = app.img.getImageCube(app.expt.filetype);
    
    % Update filename before checking the error
    app.DirectoryListBox.Items = {app.img.filename};
    app.DirectoryListBox.ItemsData = 1:length({app.img.filename});
    
    if ~isempty(ME)
        app.msgbox.Text = ME;
        return;
    end
    
else
    % Update filename
    app.DirectoryListBox.Items = {app.img.filename};
    app.DirectoryListBox.ItemsData = 1:length({app.img.filename});
end


% Update filename above listbox
app.DirectoryLabel.Text = app.img(n).label;

% If no Basis spectra selected for unmixing, do nothing
if isempty([app.spec.idx])
    app.msgbox.Text = app.msg.loadBasis;
    return;
end

% Extra check in case cube was not found
if isempty(app.img(n).cube)
    % Image wasn't loaded
    app.msgbox.Text = 'Image cube not found. This file was not loaded properly.';
    return;
end

% Analyze and unmix image if needed
if ~app.img(n).isUnmixed
    
    % Inform user
    app.msgbox.Text = app.msg.unmixing; drawnow;
    
    % Unmix image
    app.img(n) = app.img.preProcessImage(app);
    
    app.img(n) = app.img.getAmatrix(app);
    
    [app.img(n), ME] = app.img.unmixImage(app.cfg);
    if ~isempty(ME)
        app.msgbox.Text = app.msg.nogpu;
        return;
    end
    
    app.img(n) = app.img.postProcessImage(app);
    
    % Get ROI stats
    app.img(n).roiStats = getRegionStats(app.img(n), app.img(n).roi.mask);
end

app.msgbox.Text = 'Finalizing Images... Please Wait.';

% Display Evaluation time
if app.img(n).evalTime >= 1000 % If greater than 1000 ms, convert to sec
    app.UnmixingTimeEditField.Value = app.img(n).evalTime/1000;
    app.UnmixingTimeUnitLabel.Text = 's';
    
else % keep in milliseconds
    app.UnmixingTimeEditField.Value = app.img(n).evalTime;
    app.UnmixingTimeUnitLabel.Text = 'ms';
end

% Current ROI index
if app.GlobalROIMenu.Checked
    roiIdx = app.cfg.roiIdx;
else
    roiIdx = n;
end


% Display selcted image on axes1
if isempty(varargin) || any(strcmp(varargin, 'ax1'))
    
    cla(app.UIAxes1);
    switch app.SelectDisplayImageDropDown.Value
        case 1 % Basis spectra composite image
            
            specfilt = logical([app.spec.filt]); % In case filt was not a bool
            im = sum(cat(4, app.img(n).colorims.basis(specfilt).im), 4) ...
                * app.BrightnessEditField.Value; % / sum(specfilt);
            
            % Display composite image
            imshow(im, 'Parent', app.UIAxes1);
            axis(app.UIAxes1, 'tight'); 
            colorbar(app.UIAxes1, 'off');
            
        case 2 % Raw image
            
            % Raw image is normalized so don't multiply by brightness
            wlfilt = [app.img(n).filt];
            im = sum(cat(4, app.img(n).colorims.cube(wlfilt).im), 4) / sum(wlfilt);
            
            % Adjust contrast
            im = imadjust(im, stretchlim(im, 0), []);
            
            % Display raw image
            imshow(im, 'Parent', app.UIAxes1);
            axis(app.UIAxes1, 'tight'); 
            colorbar(app.UIAxes1, 'off');
            
        case 3 % Chi-squared map
            
            % Display Red. Chi-Squared Map
            imagesc(app.img(n).redchisq, 'hittest', 'off', 'Parent', app.UIAxes1);
            axis(app.UIAxes1, 'tight'); 
            colorbar(app.UIAxes1);
            app.UIAxes1.XTickLabel = [];
            app.UIAxes1.YTickLabel = [];
    end
    
    hold(app.UIAxes1, 'on');
    
    % Display ROIs
    for a = 1:length(app.img(roiIdx).roi)
        % roi = app.img(roiIdx).roi(a).mask;
        
        % Plot roi
        if isempty(app.img(roiIdx).roi(a).pts)
            % Empty roi, skip plot
        else
            switch app.img(roiIdx).roi(a).type
                case 'boundary'
                    plot(app.UIAxes1, app.img(roiIdx).roi(a).pts(:,1), app.img(roiIdx).roi(a).pts(:,2), ...
                        '.--', 'Color', app.cfg.roiColors{a}, 'MarkerSize', 15, 'LineWidth', 2);
                case 'points'
                    plot(app.UIAxes1, app.img(roiIdx).roi(a).pts(:,1), app.img(roiIdx).roi(a).pts(:,2), ...
                        '.', 'Color', app.cfg.roiColors{a}, 'MarkerSize', 1);
            end
        end
    end
    
end

% Display the selected basis image(s) on axes2
if isempty(varargin) || any(strcmp(varargin, 'ax2'))

    cla(app.UIAxes2);
    im = sum(cat(4, app.img(n).colorims.basis(m).im), 4) * app.BrightnessEditField.Value;
    imshow(im, 'Parent', app.UIAxes2);
    axis(app.UIAxes2, 'tight'); 
    colorbar(app.UIAxes2, 'off');
    
end


if isempty(varargin) || any(strcmp(varargin, 'ax3'))
    
    % Initialize vars needed for multiple axes3 options
    A  = app.img(n).A;
    wl = app.img(n).wl;
    
    % Set white color to grey for plotting
    slabels = {app.spec.label};
    scolors = {app.spec.color};
    if any(strcmp(scolors, 'White'))
        % Turn any 'White' spectra to Grey for visibility
        scolors{strcmp(scolors, 'White')} = {'Grey'};
    end
    
    % Display axes3
    cla(app.UIAxes3); 
    hold(app.UIAxes3, 'on');
    app.UIAxes3.YScale = 'linear';
    
    switch app.ChooseAnalysisDropDown.Value
        
        case 1 % Spectral Fit
            
            % Intialize data
            savg = app.img(n).roiStats(end).savg;
            sstd = app.img(n).roiStats(end).sstd;
            xavg = app.img(n).roiStats(end).xavg;
            xstd = app.img(n).roiStats(end).xstd;
            
            % Convert x to percentages
            xavg_norm = 100 * xavg / sum(xavg);
            xstd_norm = 100 * xstd;
            
            % Avg. Spectrum in ROI
            if strcmp(app.ShowErrorBarsMenu.Checked, 'on')
                % Plot error bars, but only if averaged spectra
                % Point ROI's do not have error bars...
                if isempty(sstd)
                    plot(app.UIAxes3, wl, savg, 'k');
                else
                    errorbar(app.UIAxes3, wl, savg, sstd, 'k');
                end
            else
                % Regular ol' plot
                plot(app.UIAxes3, wl, savg, 'k');
            end
            legendNames{1} = 'ROI Avg. Spec.';
            
            % Fit A * x
            plot(app.UIAxes3, wl, A*xavg, 'k--');
            legendNames{2} = 'Fit A*x';
            
            % Linewidth of each spectrum (the selected spectra get bolded)
            linewidth    = ones(app.cfg.num_spec, 1);
            linewidth(m) = 2;
            
            % Made 2 plots already, count the rest
            j = 3;
            for i = 1:app.cfg.num_spec
                if xavg(i) > 0 % Only plot if spectra is non-zero
                    
                    % Plot idv. spectra weight by unmixing coef.
                    plot(app.UIAxes3, wl, xavg(i)*A(:,i), 'Color', rgb(scolors{i}),...
                        'LineWidth', linewidth(i));
                    
                    % Keep track of legend
                    if true % isempty(xstd_norm)
                        legendNames{j} = sprintf('%s - %3.2f%%', ...
                            slabels{i}, xavg_norm(i));
                    else
                        legendNames{j} = sprintf('%s - %3.2f +/- %3.2f%%', ...
                            slabels{i}, xavg_norm(i), xstd_norm(i));
                    end
                    j = j + 1;
                end
            end
            
            app.UIAxes3.Subtitle.String = ['ROI Composition \chi^2_{avg}', ...
                    sprintf('= %4.3e', app.img(n).roiStats(end).mchisq)];
            app.UIAxes3.XLabel.String = 'Wavelength (nm)';
            app.UIAxes3.YLabel.String = 'Intensity (a.u.)';
            app.UIAxes3.XLim = [min(wl) max(wl)];
            app.UIAxes3.YLim = [0 inf];
            
            if strcmp(app.ShowLegendMenu.Checked, 'on')
                legend(app.UIAxes3, legendNames, 'Location', 'best');
                legend('boxon');
            else
                legend(app.UIAxes3, 'off');
            end
            
        case 2 % Histogram
            
            % In case multiple spectra are selected
            roi = repmat(app.img(n).roi.mask, [1 1 length(m)]);
            
            % Mask data
            data = app.img(n).x(:, :, m) .* roi;
            
            % Get rid of zeros
            data( data == 0 ) = [];
            
            % Display
            histogram(app.UIAxes3, data);
            if numel(m) == 1
                t = sprintf('Histogram of Basis Map - %s', slabels{m});
            else
                t = sprintf('Histogram of %d Basis Maps', numel(m));
            end
            app.UIAxes3.Subtitle.String = t;
            app.UIAxes3.XLabel.String = 'Intensity - 10-bit scale';
            app.UIAxes3.YLabel.String ='Counts';
            app.UIAxes3.XLim = [0 1];
            app.UIAxes3.YLim =[0 inf];
            legend(app.UIAxes3, 'off');
            
        case {3, 4} % Phasor plots
            
            % Reshape image and mask into columns of pixels
            pixels = permute(app.img(n).disp, [3 2 1]);
            %         roi    = permute(roi, [3 2 1]);
            pixels = pixels(:,:);
            %         roi    = roi(:,:);
            
            % Mask 2D image
            %         pixels = pixels(:, logical(roi));
            
            if app.ChooseAnalysisDropDown.Value == 3 % Fourier
                
                p     = app.cfg.phasor.f;
                T     = fftCoefs(pixels, p:(p+1));
                Tspec = fftCoefs(A,      p:(p+1));
                
                xlab = 'G'; ylab = 'S';
                titlestr = sprintf('Fourier Phasor plot, n = %d', app.cfg.phasor.f);
                
            elseif app.ChooseAnalysisDropDown.Value == 4 % Chebyshev
                
                p     = app.cfg.phasor.c;
                T     = fctCoefs(pixels, p:(p+1));
                Tspec = fctCoefs(A,      p:(p+1));
                
                xlab = 'U_{n}'; ylab = 'U_{n+1}';
                titlestr = sprintf('Chebyshev Phasor plot, n = %d', app.cfg.phasor.c);
                
            end
            
            % Calculate heatmap colorscheme for scatter plot
            heatmap = heatscatter(T(1,:)', T(2,:)');
            
            % Scatter plot with heatmap
            app.scat_pts = scatter(app.UIAxes3, T(1,:), T(2,:), 10, heatmap, '.'); hold on;
            
            % Plot pixels in roi
            for a = 1:length(app.img(roiIdx).roi)
                roi = app.img(roiIdx).roi(a).mask;
                roi = permute(roi, [3 2 1]);
                roi = roi(:,:);
                
                if sum(roi) < numel(roi)
                    app.scat_basis = scatter(app.UIAxes3, T(1, logical(roi)), T(2, logical(roi)), 10, app.cfg.roiColors{a}, '.');
                end
                
            end
            % Add basis spectra
            app.scat_basis = scatter(app.UIAxes3, Tspec(1,:), Tspec(2,:), 100, rgb(scolors), 'filled');
            app.UIAxes3.XLimMode = 'auto';
            app.UIAxes3.YLimMode = 'auto';
            app.UIAxes3.XLabel.String = xlab;
            app.UIAxes3.YLabel.String = ylab;
            app.UIAxes3.Subtitle.String = titlestr;
            
            if strcmp(app.ShowLegendMenu.Checked, 'on')
                legend(app.UIAxes3, {'Phasor Points', 'Basis Spectra'}, 'Location', 'best');
                colorbar(app.UIAxes3, 'on');
            else
                legend(app.UIAxes3, 'off');
                colorbar(app.UIAxes3, 'off');
            end        
    end
end

% TODO; move this code somewhere else
% Format spec roi for prism
if false
    % Intialize data
    savg = app.img(n).roiStats.savg;
    sstd = app.img(n).roiStats.sstd;
    xavg = app.img(n).roiStats.xavg;
    xstd = app.img(n).roiStats.xstd;
    npx  = app.img(n).roiStats.npx;
    npxvec = repmat(npx, [1 length(wl)]);
    j = 1;
%     assignin('base', 'A', A);
%     assignin('base', 'xavg', xavg);
%     assignin('base', 'savg', savg);
    
    for i = 1:3:((app.cfg.num_spec+2)*3)
        if i == 1
            
            specdata2prism(:, i)     = savg;
            specdata2prism(:, i + 1) = sstd / sqrt(npx);
            specdata2prism(:, i + 2) = npxvec;
            
        elseif i == 4
            
            specdata2prism(:, i)     = A*xavg;
            specdata2prism(:, i + 1) = A*xstd;
            specdata2prism(:, i + 2) = npxvec;
            
        else
            
            specdata2prism(:, i)     = xavg(j)*A(:,j);
            specdata2prism(:, i + 1) = xstd(j)*A(:,j);
            specdata2prism(:, i + 2) = npxvec;

            j = j + 1;
        end
        
        
    end
    
    % Sum-to-one norm
    xsum = sum(xavg);
    xavg = 100*xavg/xsum;
    xstd = 100*xstd/xsum;
    xdata2prism = zeros(1, length(xavg)*3);
    xdata2prism(1:3:end) = xavg;
    xdata2prism(2:3:end) = xstd;
    xdata2prism(3:3:end) = npx;
    
% %     xdata2prism = [xavg(:) xstd(:) repmat(npx, [cfg.num_spec 1])];
    
    % Pass data to workspace
    assignin('base', 'specdata2prism', specdata2prism);
    assignin('base', 'xdata2prism',    xdata2prism);
%     assignin('base', 'A', A);
%     assignin('base', 'xavg', xavg);
%     assignin('base', 'savg', savg);
    
end


app.msgbox.Text = 'Ready.';
drawnow;