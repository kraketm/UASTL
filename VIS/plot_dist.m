function plot_dist(data, type, opts)
% This function PLOT_DIST creates the individual plot of the distribution 
% for given mean and covariance matrix.
%
%
%  [] = PLOT_DIST(data, type, opts)
%
% Input:
%   * data: data struct with data.mu and data.Sigma
%   * type: one of "comb", "isoband", or "spaghetti".
%           Here, "comb" combines isobands and spaghetti samples.
%   * opts: optional input (to individualize plotting and exporting)
%       - lineWidth {mustBeNumeric} = 1
%       - samplesColored (1,1) logical = false
%       - nmbsamples (1,1) {mustBeNumeric} = 5 - Number of samples 
%       - smplWidth (1,1) = 0.5
%       - ylim (1,2) {mustBeNumeric} = [0 0] - scope of yaxis
%       - pbaspect (3,1) {mustBeNumeric} = [6 1 1] - pbaspect of axis plot
%       - export (1,:) struct = struct - can contain the fields 
%       - export_path, exportName, export, pbaspect, yaxis, dashedLines, and format
%
% Output:
%   * plots into current figure handle, if no figure is open, new figure 


arguments
    data (1,1) struct
    type {mustBeMember(type,["comb","isoband","spaghetti"])} = "comb"
    opts.lineWidth {mustBeNumeric} = 1;
    opts.samplesColored (1,1) logical = false;
    opts.nmbsamples (1,1) {mustBeNumeric} = 5;
    opts.smplWidth (1,1) = 0.5;
    opts.ylim (1,2) {mustBeNumeric} = [0,0];
    opts.pbaspect (3,1) {mustBeNumeric} = [6 1 1];
    opts.export (1,:) struct = struct;
end
lineFact = opts.lineWidth;

if isempty(fieldnames(opts.export))
    opts.export.export = false;
else
    mustBeMember(fieldnames(opts.export), ["exportPath", "exportName", "export", ...
        "pbaspect", "yaxis", "dashedLines", "format"])
    if ~isfield(opts.export, "exportPath")
        opts.export.exportPath = "";
    end
    if ~isfield(opts.export, "exportName")
        opts.export.exportName = "plot_dist";
    end
    if isfield(opts.export, "pbaspect")
        opts.pbaspect = opts.export.pbaspect;
    end
end

if ~(isfield(data, "mu")); error 'mu not defined.'; end

if ~isfield(data, "Sigma") && ~isfield(data, "sigma")
    error 'Sigma and sigma not defined.'
end

if ~isfield(data, "sigma"); data.sigma = real(sqrt(diag(data.Sigma))); end

if ~isfield(data, "samples")
    data.samples = mvnrnd(data.mu,data.Sigma,opts.nmbsamples)';
end

% Sigma levels: 0%, 50%, 99%
sigmalvl = [0, 0.674490, 2.575829];
set(gcf,'Renderer','opengl');

% Define own predefined colors
OwnColors = [228,026,028;055,126,184;077,175,074;152,078,163;...
    255,127,000;255,255,051;166,086,040;251,180,174;179,205,227;...
    204,235,197;222,203,228;254,217,166;255,255,204;229,216,189]./255;
col(:,1) = linspace(OwnColors(2,1), 1, numel(sigmalvl)+2);
col(:,2) = linspace(OwnColors(2,2), 1, numel(sigmalvl)+2);
col(:,3) = linspace(OwnColors(2,3), 1, numel(sigmalvl)+2);

if isequal(opts.ylim,[0,0])
    axbounds = [1, length(data.mu), min(data.mu - ceil(max(sigmalvl)+1)*data.sigma), max(data.mu + ceil(max(sigmalvl)+1)*data.sigma)];
else
    axbounds = [1, length(data.mu), opts.ylim(1), opts.ylim(2)];
end

if type == "isoband"
    hold on
    x=1:length(data.mu);
    for j=length(sigmalvl):-1:2
        posCont = data.mu + sigmalvl(j)*data.sigma;
        negCont = data.mu - sigmalvl(j)*data.sigma;
        for i = 1:length(data.mu)-1
            xp = [x(i) x(i+1) x(i+1) x(i)];
            yp = [negCont(i) negCont(i+1) posCont(i+1) posCont(i)];
            patch(xp, yp, col(j,:), 'EdgeColor','none', 'FaceAlpha',.5);
        end
    end
    plot(x,data.mu, 'Color', col(1,:), 'LineWidth',lineFact*1.5);
    hold off
elseif type == "spaghetti"
    initialColorOrder = [...
        0.9290,    0.6940,    0.1250; 0.4660,    0.6740,    0.1880; ...
        0.4940,    0.1840,    0.5560; 0.6350,    0.0780,    0.1840; ...
        0.8500,    0.3250,    0.0980; 0 ,   0.4470 ,   0.7410; ...
        0.3010,    0.7450,    0.9330];
    if opts.nmbsamples > length(initialColorOrder)
        initialColorOrder = repmat(initialColorOrder, ...
            [floor(opts.nmbsamples/length(initialColorOrder)),1]);
        initialColorOrder = [initialColorOrder; ...
            initialColorOrder(1:mod(opts.nmbsamples, length(initialColorOrder)),:)];
    else
        initialColorOrder = [initialColorOrder(1:opts.nmbsamples,:)];
    end
    colors = initialColorOrder;
    if size(data.samples,2) == 2* opts.nmbsamples
        hold on
        samples2plot = data.samples(:,1:opts.nmbsamples);
        samples2plot_shift = data.samples(:,opts.nmbsamples+1:2*opts.nmbsamples);

        if ~opts.samplesColored
            h2a = plot(samples2plot, 'Color',OwnColors(2,:), 'LineWidth',lineFact*opts.smplWidth); 
            h2a_shift = plot(samples2plot_shift, 'Color',OwnColors(2,:), 'LineWidth',lineFact*opts.smplWidth*0.2, 'LineStyle','-'); 
            h2a_shift_dot = plot(samples2plot_shift, 'Color',OwnColors(2,:), 'LineWidth',lineFact*opts.smplWidth, 'LineStyle',':'); 
        else
            for i=1:opts.nmbsamples
                plot(samples2plot(:,i), 'LineWidth',lineFact*opts.smplWidth, 'Color',colors(i,:)); 
                plot(samples2plot_shift(:,i), 'LineWidth',lineFact*opts.smplWidth*0.2, 'Color',colors(i,:), 'LineStyle','-'); 
                plot(samples2plot_shift(:,i), 'LineWidth',lineFact*opts.smplWidth, 'Color',colors(i,:), 'LineStyle',':'); 
            end
        end
        if ~opts.samplesColored && opts.nmbsamples>1
            for i=1:numel(h2a)
                h2a(i).Color(4) = 0.5;
                h2a_shift(i).Color(4) = 0.5;
                h2a_shift_dot(i).Color(4) = 0.5;
            end
        end
    else
        hold on
        samples2plot = data.samples(:,1:opts.nmbsamples);

        if ~opts.samplesColored
            h2a = plot(samples2plot, 'Color',OwnColors(2,:), 'LineWidth',lineFact*opts.smplWidth); 
        else
            for i=1:opts.nmbsamples
                plot(samples2plot(:,i), 'LineWidth',lineFact*opts.smplWidth, 'Color',colors(i,:)); 
            end
        end
        if ~opts.samplesColored && opts.nmbsamples>1
            for i=1:numel(h2a)
                h2a(i).Color(4) = 0.5;
            end
        end
    end
elseif type == "comb"
    plot_dist(data, "isoband", samplesColored = opts.samplesColored, ylim=opts.ylim)
    
    hold on
    if isfield(data, "sigmaNew")
        %sigmalvl(2)
        plot(1:length(data.mu),data.mu+sigmalvl(2)*data.sigmaNew, 'Color', [col(2,:),0.5], 'LineWidth',lineFact*0.25, 'LineStyle','-');
        plot(1:length(data.mu),data.mu-sigmalvl(2)*data.sigmaNew, 'Color', [col(2,:),0.5], 'LineWidth',lineFact*0.25, 'LineStyle','-');
        %sigmalvl(3)
        plot(1:length(data.mu),data.mu+sigmalvl(3)*data.sigmaNew, 'Color', [col(3,:),0.5], 'LineWidth',lineFact*0.25, 'LineStyle','-');
        plot(1:length(data.mu),data.mu-sigmalvl(3)*data.sigmaNew, 'Color', [col(3,:),0.5], 'LineWidth',lineFact*0.25, 'LineStyle','-');
        %sigmalvl(2)
        plot(1:length(data.mu),data.mu+sigmalvl(2)*data.sigmaNew, 'Color', [col(2,:),0.5], 'LineWidth',lineFact*0.5, 'LineStyle',':');
        plot(1:length(data.mu),data.mu-sigmalvl(2)*data.sigmaNew, 'Color', [col(2,:),0.5], 'LineWidth',lineFact*0.5, 'LineStyle',':');
        %sigmalvl(3)
        plot(1:length(data.mu),data.mu+sigmalvl(3)*data.sigmaNew, 'Color', [col(3,:),0.5], 'LineWidth',lineFact*0.5, 'LineStyle',':');
        plot(1:length(data.mu),data.mu-sigmalvl(3)*data.sigmaNew, 'Color', [col(3,:),0.5], 'LineWidth',lineFact*0.5, 'LineStyle',':');
    end
    plot(1:length(data.mu),data.mu, 'Color', col(1,:), 'LineWidth',lineFact*1);
    plot_dist(data, "spaghetti", samplesColored = opts.samplesColored, ylim=opts.ylim, nmbsamples = opts.nmbsamples, lineWidth=lineFact*1, smplWidth = opts.smplWidth)
    if opts.export.export
        opts.axbounds = axbounds;
        exportSubPlot(opts)
    end
    return;
end

if opts.export.export
    opts.axbounds = axbounds;
    exportSubPlot(opts)
end
if axbounds(3)~=axbounds(4)
    axis(axbounds);
end
end % Function plot_dist


function exportSubPlot(opts)
% This function EXPORTSUBPLOT uses the current axis handle and exports it
% using the specified export options (opts.export.*).
if ~isfield(opts.export, "format"); opts.export.format = "SUBs"; end
ax = gca;
set(ax,'XTick',[], 'XColor', [1 1 1])
set(ax,'YTick',[], 'YColor', [1 1 1])
set(ax, xlim=opts.axbounds(1:2))
if any(opts.axbounds(3:4)~=[0 0])
set(ax, ylim=opts.axbounds(3:4))
end
set(ax, 'Color', 'none')
set(ax,'Visible','on')
pbaspect(ax,opts.pbaspect)
%box on
colorbar off
if opts.export.format == "SUBs"
    if opts.export.export
        exportgraphics(ax, strcat(opts.export.exportPath,opts.export.exportName,".pdf"),...
            'ContentType','vector', 'BackgroundColor','none', 'Resolution',900);
    end
elseif opts.export.format == "FIG"
    set(ax,'XColor', [0 0 0])
    set(ax, 'YColor', [0 0 0])
    box on
    if isfield(opts.export, "dashedLines")
        set(ax, 'YTick',[opts.export.dashedLines])
    end
end
end % Function exportSubPlot