function displayImage(app)
%% displayImage - update the GUI display

% If no experiment, do nothing
if ~app.expt.loaded
    app.msgbox.text = app.msg.loadexpt;
    return;
end

% Current selected image
n = app.state.img;

if strcmp(app.img(n).filename, 'NOT_LOADED')
    % File was not loaded properly the first time, skip it
    app.msgbox.text = 'File not loaded properly.';
    return;
end

% Local varaibles
cfg   = app.cfg;
state = app.state;
spec  = app.spec;
expt  = app.expt;
m     = state.spec;

% Check if image has been loaded
if ~app.img(n).isLoaded
    % If not, load image
    
    % Inform user
    app.msgbox.text = app.msg.loadingcube;
    
    [app.img(n), ME] = getImageCube(app.img(n), expt.filetype);
    
    % Update filename before checking the error
    app.img_listbox.String = {app.img.filename};
    
    if ~isempty(ME)
        app.msgbox.text = ME;
        return;
    end
    
else
    % Update filename
    app.img_listbox.String = {app.img.filename}; 
end



% Update filename above listbox
app.text_filename.String = app.img(n).label;

% If no Basis spectra selected for unmixing, do nothing
if isempty([spec.idx])
    app.msgbox.text = app.msg.loadBasis;
    return;
end

% Extra check incase cube was not found
if isempty(app.img(n).cube)
    % Image wasn't loaded
    app.msgbox.text = 'This file was not loaded properly.';
    return;
end

% Analyze and unmix image if needed
if ~app.img(n).isUnmixed
    
    % Inform user
    app.msgbox.text = app.msg.unmixing; drawnow;
    
    % Unmix image
    [app.img(n), ME] = unmixImage(app.img(n), app);
    if ~isempty(ME)
        app.msgbox.text = app.msg.nogpu;
        return;
    end
    
    % Get ROI stats
    app.img(n).roiStats = getRegionStats(app.img(n), app.img(n).roi.mask);
end

app.msgbox.text = 'Finalizing Images... Please Wait.';

% Format evaluation time
if app.img(n).evalTime >= 1000 % If greater than 1000 ms, convert to sec
    timestr = sprintf('%3.1f s', app.img(n).evalTime/1000);
else % keep in milliseconds
    timestr = sprintf('%3.1f ms', app.img(n).evalTime);
end

% Display time
app.text_time.String = timestr;

% Display selcted image on axes1
axes(app.axes1); cla;
switch app.display_ax1_popupmenu.Value
    case 1 % Basis spectra composite image
        
        specfilt = logical([app.spec.filt]); % In case filt was not a bool
        im = sum(cat(4, app.img(n).colorims.basis(specfilt).im), 4) ...
             * state.brightness; % / sum(specfilt);

        % Display composite image
        imObj = imshow(im);
        axis tight; colorbar('off');

    case 2 % Raw image
        
        % Raw image is normalized so don't multiply by brightness
        wlfilt = [app.img(n).filt];
        im = sum(cat(4, app.img(n).colorims.cube(wlfilt).im), 4) / sum(wlfilt);
        
        % Adjust contrast
        im = imadjust(im, stretchlim(im, 0), []);
        
        % Display raw image
        imObj = imshow(im);
        axis tight; colorbar('off');
        
    case 3 % Chi-squared map
                
        % Display Red. Chi-Squared Map
        imObj = imagesc(app.img(n).redchisq, 'hittest', 'off');
        axis tight;
        app.axes1.XTickLabel = [];
        app.axes1.YTickLabel = [];
        colorbar;
end

% Give the image context menu properties
set(imObj, 'uicontextmenu', app.axes1_cmenu);
% app.imObj = imObj;
hold on;

% Current ROI index
if state.globalroi
    roiIdx = cfg.roiIdx;
else
    roiIdx = n;
end

% Display ROIs
for a = 1:length(app.img(roiIdx).roi)
    roi = app.img(roiIdx).roi(a).mask;
    
    % Plot roi
    if isempty(app.img(roiIdx).roi(a).pts)
        % Empty roi, skip plot
    else
        switch app.img(roiIdx).roi(a).type
            case 'boundary'
                plot(app.img(roiIdx).roi(a).pts(:,1), app.img(roiIdx).roi(a).pts(:,2), ...
                    '.--', 'Color', cfg.roiColors{a}, 'MarkerSize', 15, 'LineWidth', 2);
            case 'points'
                plot(app.img(roiIdx).roi(a).pts(:,1), app.img(roiIdx).roi(a).pts(:,2), ...
                    '.', 'Color', cfg.roiColors{a}, 'MarkerSize', 1);
        end
    end
end
 
% Display the selected basis image(s) on axes2
axes(app.axes2); cla;
im = sum(cat(4, app.img(n).colorims.basis(m).im), 4) * state.brightness;
imObj = imshow(im);
set(imObj, 'uicontextmenu', app.axes1_cmenu);

% Initialize vars needed for multiple axes3 options
A  = app.img(n).A;
wl = app.img(n).wl;

% Set white color to grey for plotting
slabels = {spec.label};
scolors = {spec.color};
try
    scolors{strcmp(scolors, 'White')} = 'Grey';
catch
end

% Display axes3
axes(app.axes3); cla; hold on;
set(gca, 'YScale', 'lin');

switch app.display_ax3_popupmenu.Value
    
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
        if state.errorbars
            % Plot error bars, but only if averaged spectra
            % Point ROI's do not have error bars...
            if isempty(sstd)
                plot(wl, savg, 'k');
            else
                errorbar(wl, savg, sstd, 'k');
            end
        else
            % Regular ol' plot
            plot(wl, savg, 'k');
        end
        legendNames{1} = 'ROI Avg. Spec.';
        
        % Fit A * x
        plot(wl, A*xavg, 'k--');
        legendNames{2} = 'Fit A*x';
        
        % Linewidth of each spectrum (the selected spectra get bolded)
        linewidth    = ones(cfg.num_spec, 1);
        linewidth(m) = 2;
        
        % Made 2 plots already, count the rest
        j = 3;
        for i = 1:cfg.num_spec
            if xavg(i) > 0 % Only plot if spectra is non-zero
                
                % Plot idv. spectra weight by unmixing coef.
                plot(wl, xavg(i)*A(:,i), 'Color', rgb(scolors{i}),...
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
        
        t = ['ROI Composition \chi^2_{avg}' sprintf('= %4.3e', app.img(n).roiStats(end).mchisq)];
        title(t);
        xlabel('Wavelength (nm)');
        ylabel('Intensity (a.u.)');
        xlim([min(wl) max(wl)]);
        ylim([0 inf]);
        
        if state.legend
            legend(legendNames, 'Location', 'best');
            legend('boxon');
        else
            legend off
        end
        
    case 2 % Histogram
        
        % In case multiple spectra are selected
        roi = repmat(roi, [1 1 length(m)]);
        
        % Mask data
        data = app.img(n).x(:, :, m) .* roi;
        
        % Get rid of zeros
        data( data == 0 ) = [];
        
        % Display
        histogram(data);
        if numel(m) == 1
            t = sprintf('Histogram of Basis Map - %s', slabels{m});
        else
            t = sprintf('Histogram of %d Basis Maps', numel(m));
        end
        title(t);
        xlabel('Intensity - 10-bit scale');
        ylabel('Counts');
        xlim([0 1]);
        ylim([0 inf]);
        
    case {3, 4} % Phasor plots
        
        % Reshape image and mask into columns of pixels
        pixels = permute(app.img(n).disp, [3 2 1]);
%         roi    = permute(roi, [3 2 1]);
        pixels = pixels(:,:);
%         roi    = roi(:,:);
        
        % Mask 2D image
%         pixels = pixels(:, logical(roi));
                
        if app.display_ax3_popupmenu.Value == 3 % Fourier
            
            p = cfg.phasor.f;
            T     = fftCoefs(pixels, p:(p+1));
            Tspec = fftCoefs(A,      p:(p+1));
            
            xlab = 'G'; ylab = 'S'; 
            titlestr = sprintf('Fourier Phasor plot, n = %d', cfg.phasor.f);
              
        elseif app.display_ax3_popupmenu.Value == 4 % Chebyshev
            
            p     = cfg.phasor.c;
            T     = fctCoefs(pixels, p:(p+1));
            Tspec = fctCoefs(A,      p:(p+1));
            
            xlab = 'U_{n}'; ylab = 'U_{n+1}';
            titlestr = sprintf('Chebyshev Phasor plot, n = %d', cfg.phasor.c);
            
        end
        
        % Calculate heatmap colorscheme for scatter plot
        heatmap = heatscatter(T(1,:)', T(2,:)');
        
        % Scatter plot with heatmap
        app.scat_pts = scatter(T(1,:), T(2,:), 10, heatmap, '.'); hold on;
        
        % Plot pixels in roi
        for a = 1:length(app.img(roiIdx).roi)
            roi = app.img(roiIdx).roi(a).mask;
            roi = permute(roi, [3 2 1]);
            roi = roi(:,:);
        
            if sum(roi) < numel(roi)
                app.scat_basis = scatter(T(1, logical(roi)), T(2, logical(roi)), 10, cfg.roiColors{a}, '.');
            end
        
        end
        % Add basis spectra
        app.scat_basis = scatter(Tspec(1,:), Tspec(2,:), 100, rgb(scolors), 'filled');
        set(gca, 'XLimMode', 'auto'); set(gca, 'YLimMode', 'auto');
%         xlim([0 inf]); ylim([0 inf]);
        xlabel(xlab); ylabel(ylab); title(titlestr);
        
        if state.legend
            legend({'Phasor Points', 'Basis Spectra'}, 'Location', 'best');
            colorbar;
        else
            legend off
            colorbar off;
        end
 
end

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
    
    for i = 1:3:((cfg.num_spec+2)*3)
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


app.msgbox.text = 'Ready.';
drawnow;