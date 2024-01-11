function yLTSTR = plot_distributionmtx(yLTSTR, numPeriods, type, opts)
% This function PLOT_DISTRIBUTIONMTX creates the main visualization plot
% of the distribution for given mean and covariance matrix.
% It is highly dependent from the function plot_dist.
%
% Plot distribution for time series depending the Gaussian distributed 
% variable yLTSTR = [y, LT, ST_1, ..., ST_L, R] with
%   * yLTSTR.mu: mean vector 
%   * yLTSTR.Sigma: covariance matrix
%
% The theory is described in the related manuscript "Uncertainty-Aware
% Seasonal-Trend Decomposition Based on Loess".
%
%
% yLTSTR = PLOT_DISTRIBUTIONMTX(yLTSTR, numPeriods, type, opts)
%
% Input:
%   * yLTSTR: struct with mean and Sigma containing the data y, trend LT, 
%             seasonal ST, and residual R components computed by UASTL
%       - yLTSTR.mu: mean vector with components
%       - yLTSTR.Sigma: covariance matrix with components
%   * numPeriods: Number of periods to subdevide the yLTSTR cell data
%   * type: one of "comb", "isoband", or "spaghetti"
%           "comb" plots isobands and spaghetti samples
%   * opts: optional input with additional VIS techniques
%       - timeStamp (:,1) double = []
%       - nmbsamples (1,1) {mustBeInteger} = 5
%       - samplesColored (1,1) logical = false
%       - lineWidth {mustBeNumeric} = 1;
%       - plotCov (1,1) logical = false
%       - plotCor (1,1) logical = false
%       - plotCorLength (1,:) struct = struct
%       - coPoint (1,1) {mustBeInteger} = 0
%       - delta (:,1) double = []
%       - export
%
% Output: builds Figure(s) containing the desired VIS techniques for UASTL
%   * yLTSTR: Same yLTSTR but with added field CorMat if correlation matrix
%             is computed (if coPoint or plotCor options are set)
%

arguments
    yLTSTR (1,:) struct
    numPeriods (1,1) {mustBeInteger}
    type {mustBeMember(type,["comb","isoband", "spaghetti"])}
    opts.timeStamp (:,1) double = []
    opts.nmbsamples (1,1) {mustBeInteger} = 5
    opts.samplesColored (1,1) logical = false
    opts.lineWidth {mustBeNumeric} = 1;
    opts.discrNmb (1,1) {mustBeInteger} = 13;
    opts.plotCov (1,1) logical = false
    opts.plotCor (1,1) logical = false
    opts.plotCorLength (1,1) logical = false
    opts.coPoint (1,1) {mustBeInteger} = 0
    opts.delta (:,1) double = []
    opts.export (1,:) struct = struct
end

fprintf('############################### Visualization ################################### \n')
fprintf('##### Uncertainty-Aware Seasonal-Trend Decomposition Based on Loess (UASTL) ##### \n')
fprintf('################################################################################# \n')
fprintf('##### \n')
fprintf('##### Number of Periods:\t numPeriods = %i \n',numPeriods)
fprintf('##### \n')
fprintf('##### Plotting Type:\t \t type = %s \n',string(type))
fprintf('##### \t ├── Number of Spaghetthi: \t nmbsamples = %i \n',opts.nmbsamples)
fprintf('##### \t ├── Spaghetthi colored:\t samplesColored = %s \n',string(opts.samplesColored))
fprintf('##### \t ├── Set LineWidth Factor:\t lineWidth = %.1d \n',opts.lineWidth)
fprintf('##### \t └── Discretization Steps:\t discrNmb = %i \n',opts.discrNmb)
fprintf('##### \n')
if ~isempty(opts.timeStamp)
    fprintf('##### Plotting Interval:\t timeStamp = [%i,%i] \n',opts.timeStamp(1),opts.timeStamp(2))
else
    fprintf('##### Plotting Interval:\t Full Time Domain.\n')
end
fprintf('##### \n')
fprintf('##### Visualization Techniques:\n')
fprintf('##### \t ├── Covariance Matrix:\t\t plotCov = %s \n', string(opts.plotCov))
fprintf('##### \t ├── Correlation Matrix:\t plotCor = %s \n',string(opts.plotCor))
fprintf('##### \t ├── Correlation Plot:\t \t plotCorLength = %s \n',string(opts.plotCorLength))
if opts.coPoint~=0
    fprintf('##### \t └── Local Corr. Details:\t coPoint = %i \n', opts.coPoint)
else
    fprintf('##### \t └── No Local Corr. Details:\t coPoint = %i \n',opts.coPoint)
end
fprintf('##### \n')
if isfield(opts.export, 'export') && (opts.export.export == true) && isfield(opts.export, 'exportPath')
    fprintf('##### Subplot Export into Folder: \t export.exportPath = %s \n', string(opts.export.exportPath))
elseif isfield(opts.export, 'export') && (opts.export.export == true)
    fprintf('##### Subplots Export into Current Folder')
else
    fprintf('##### Subplots not exported: \t\t export.export = false\n')
end
fprintf('##### \n')
fprintf('################################################################################# \n')
fprintf('##### \n')
fprintf('##### Plotting START \n')
fprintf('##### \n')

%%%% checks and setup of (optional) variables
% yLTSTR.mu & yLTSTR.Sigma
if (numel(yLTSTR.mu) == size(yLTSTR.Sigma, 1)) && (numel(yLTSTR.mu) == size(yLTSTR.Sigma, 2))
    ylen = numel(yLTSTR.mu)/(3+numPeriods);
else
    error 'Size of mu and Sigma are not concise'
end

% export
if ~isempty(fieldnames(opts.export))
    mustBeMember(fieldnames(opts.export), ["exportPath", "exportName", "export", "pbaspect", "yaxis", "dashedLines", "format"])

    if ~isfield(opts.export, "exportPath") && opts.export.export
        opts.export.exportPath = ""; 
    end
    if isfield(opts.export, "format"); mustBeMember(opts.export.format, ["FIG", "SUBs"]); end

    if isfield(opts.export, "exportName")
        exportName = opts.export.exportName;
    else
        exportName = "plot_dist";
    end
else
    opts.export.export = false;
    opts.export.format = "FIG";
    exportName = "";
end

% discrNmb
if opts.discrNmb > 257
    opts.discrNmb = 257;
    warning '"nmbColors" set too large, setting to 257'
elseif mod(opts.discrNmb,2) == 0
    opts.discrNmb = opts.discrNmb+1;
    warning '"nmbColors" set to the next higher odd number.'
end


if ~isempty(opts.timeStamp) && (opts.timeStamp(2) <= opts.timeStamp(1))
    warning 'TimeStamp not correctly assigned. Ignoring input.'
    opts.timeStamp = [];
end

if isempty(opts.timeStamp)
    TmeStmp = [1, ylen];
    opts.StmpLength = ylen;
    opts.stamp_indices = (1:length(yLTSTR.mu));
    opts.timeStamp = [1, ylen];
else
    TmeStmp = opts.timeStamp;
    opts.StmpLength = TmeStmp(2)-TmeStmp(1)+1;
    stamp_indices = opts.timeStamp(1):opts.timeStamp(2);
    stamp_indices = [stamp_indices, ylen+opts.timeStamp(1):ylen+opts.timeStamp(2)];
    for i=1:numPeriods
        stamp_indices = [stamp_indices, (i+1)*ylen+opts.timeStamp(1):(i+1)*ylen+opts.timeStamp(2)];
    end
    opts.stamp_indices = [stamp_indices, (numPeriods+2)*ylen+opts.timeStamp(1):(numPeriods+2)*ylen+opts.timeStamp(2)];
end


% Picking global samples if comb or spaghetti plot
if type == "comb" || type == "spaghetti"
    if ~isempty(opts.delta) && isfield(yLTSTR, "AHat")
        % SHIFTED samples via delta vector

        T = cholcov(yLTSTR.Sigma(1:ylen, 1:ylen))';
        sampOrig = randn(opts.nmbsamples, size(T,2))';

        AHat_EHat = yLTSTR.AHat;
        T_sampOrig = T*sampOrig;
        
        samp = AHat_EHat*T_sampOrig + yLTSTR.mu;
        sampNeu = samp + AHat_EHat*diag(opts.delta-1)*T_sampOrig;

        yLTSTR.samples = [samp,sampNeu];

        sigmaNew = sqrt( sum((AHat_EHat*diag(opts.delta)*T).^2,2) );

        y.samples = yLTSTR.samples(TmeStmp(1):TmeStmp(2),:);
        y.sigmaNew = sigmaNew(TmeStmp(1):TmeStmp(2));

        LT.samples = yLTSTR.samples(ylen+TmeStmp(1):ylen+TmeStmp(2),:);
        LT.sigmaNew = sigmaNew(ylen+TmeStmp(1):ylen+TmeStmp(2));
        
        for i=1:numPeriods
            ST.samples(:,:,i) = yLTSTR.samples((1+i)*ylen+TmeStmp(1):(1+i)*ylen+TmeStmp(2),:);
            ST.sigmaNew(:,i) = sigmaNew((1+i)*ylen+TmeStmp(1):(1+i)*ylen+TmeStmp(2));
        end
        R.samples = yLTSTR.samples((2+numPeriods)*ylen+TmeStmp(1):(2+numPeriods)*ylen+TmeStmp(2),:);
        R.sigmaNew = sigmaNew((2+numPeriods)*ylen+TmeStmp(1):(2+numPeriods)*ylen+TmeStmp(2),:);
    else
        if ~isempty(opts.delta)
            warning 'No yLTSTR.AHat given, therefore going on without delta.'
        end
        yLTSTR.samples = mvnrnd(yLTSTR.mu,yLTSTR.Sigma,opts.nmbsamples)';
        y.samples = yLTSTR.samples(TmeStmp(1):TmeStmp(2),:);
        LT.samples = yLTSTR.samples(ylen+TmeStmp(1):ylen+TmeStmp(2),:);
        for i=1:numPeriods
            ST.samples(:,:,i) = yLTSTR.samples((1+i)*ylen+TmeStmp(1):(1+i)*ylen+TmeStmp(2),:);
        end
        R.samples = yLTSTR.samples((2+numPeriods)*ylen+TmeStmp(1):(2+numPeriods)*ylen+TmeStmp(2),:);
    end
end


y.mu = yLTSTR.mu(TmeStmp(1):TmeStmp(2));
y.sigma = sqrt(max(diag(yLTSTR.Sigma(TmeStmp(1):TmeStmp(2), TmeStmp(1):TmeStmp(2))),0));
y.Sigma = yLTSTR.Sigma(TmeStmp(1):TmeStmp(2), TmeStmp(1):TmeStmp(2));

LT.mu = yLTSTR.mu(ylen+TmeStmp(1): ylen+TmeStmp(2));
LT.sigma = sqrt(max(diag(yLTSTR.Sigma(ylen+TmeStmp(1): ylen+TmeStmp(2), ylen+TmeStmp(1): ylen+TmeStmp(2))),0));
LT.Sigma = yLTSTR.Sigma(ylen+TmeStmp(1): ylen+TmeStmp(2), ylen+TmeStmp(1): ylen+TmeStmp(2));

for i=1:numPeriods
    ST.mu(:,i) = yLTSTR.mu((1+i)*ylen+TmeStmp(1): (1+i)*ylen+TmeStmp(2));
    ST.sigma(:,i) = sqrt(max(diag(yLTSTR.Sigma((1+i)*ylen+TmeStmp(1): (1+i)*ylen+TmeStmp(2), (1+i)*ylen+TmeStmp(1): (1+i)*ylen+TmeStmp(2))),0));
    ST.Sigma(:,:,i) = yLTSTR.Sigma((1+i)*ylen+TmeStmp(1): (1+i)*ylen+TmeStmp(2),(1+i)*ylen+TmeStmp(1): (1+i)*ylen+TmeStmp(2));
end

R.mu = yLTSTR.mu((2+numPeriods)*ylen+TmeStmp(1): (2+numPeriods)*ylen+TmeStmp(2));
R.sigma = sqrt(max(diag(yLTSTR.Sigma((2+numPeriods)*ylen+TmeStmp(1): (2+numPeriods)*ylen+TmeStmp(2),(2+numPeriods)*ylen+TmeStmp(1): (2+numPeriods)*ylen+TmeStmp(2))),0));
R.Sigma = yLTSTR.Sigma((2+numPeriods)*ylen+TmeStmp(1): (2+numPeriods)*ylen+TmeStmp(2),(2+numPeriods)*ylen+TmeStmp(1): (2+numPeriods)*ylen+TmeStmp(2));

% SET FIGURE
figure
% get(0, 'Screensize') -> [1 1 1710 1112]
set(gcf, 'Position', [1 1 1710 1112]);
sgtitle('Uncertainty-Aware Seasonal-Trend Decomposition')

nmbColors = opts.discrNmb;
colors = ColorLUT(1);
colors = colors(round(linspace(1,257,nmbColors)),:);

if opts.coPoint==0
    helperCoDep = 0;
elseif (abs(opts.coPoint) > length(yLTSTR.mu))
    warning 'Covariance point set too high, plotting without dependency plot.\n'
    helperCoDep = 0;
elseif any(opts.stamp_indices==opts.coPoint)
    helperCoDep = 1;
else
    helperCoDep = 1;
    warning 'Covariance point is out of timeStamp-Range.\n'
end


if helperCoDep==1
    Interactiv_Pointer = opts.coPoint;
    %%% We choose CorMat for background plot.
    yLTSTR = corComp(yLTSTR);
    plotBack = yLTSTR.corMat(Interactiv_Pointer,:);

end
% Get the ylims for each part by (mean+sigmalvlmax)/3*4 to add below
sigmalvl = [0, 0.674490, 2.575829];

part_factor = 5;
y.lims = [1, length(y.mu), min(y.mu - ceil(max(sigmalvl))*y.sigma), max(y.mu + ceil(max(sigmalvl))*y.sigma)];
y.lims(3) = y.lims(3) - helperCoDep * 1/part_factor * (y.lims(4)-y.lims(3));
maxyheight(1) = ((y.lims(4)-y.lims(3))/(part_factor + 1))/2;
ymid(1) = y.lims(3) + maxyheight(1);

LT.lims = [1, length(LT.mu), min(LT.mu - ceil(max(sigmalvl))*LT.sigma), max(LT.mu + ceil(max(sigmalvl))*LT.sigma)];
LT.lims(3) = LT.lims(3) - helperCoDep * 1/part_factor * (LT.lims(4)-LT.lims(3));
maxyheight(2) = ((LT.lims(4)-LT.lims(3))/(part_factor + 1))/2;
ymid(2) = LT.lims(3) + maxyheight(2);
for i=1:numPeriods
    ST.lims(:,i) = [1, length(ST.mu(:,i)), min(ST.mu(:,i) - ceil(max(sigmalvl))*ST.sigma(:,i)), max(ST.mu(:,i) + ceil(max(sigmalvl))*ST.sigma(:,i))];
    ST.lims(3,i) = ST.lims(3,i) - helperCoDep * 1/part_factor * (ST.lims(4,i)-ST.lims(3,i));
    maxyheight(i+2) = ((ST.lims(4,i)-ST.lims(3,i))/(part_factor + 1))/2;
    ymid(i+2) = ST.lims(3,i) + maxyheight(i+2);
end

R.lims = [1, length(R.mu), min(R.mu - ceil(max(sigmalvl))*R.sigma), max(R.mu + ceil(max(sigmalvl))*R.sigma)];
R.lims(3) = R.lims(3) - helperCoDep * 1/part_factor * (R.lims(4)-R.lims(3));
maxyheight(numPeriods+3) = ((R.lims(4)-R.lims(3))/(part_factor + 1))/2;
ymid(numPeriods+3) = R.lims(3) + maxyheight(numPeriods+3);

x=linspace(1,ylen,ylen);

% Approx. 16:9 pbaspect of plot:
pbxLength = (3+numPeriods + 1)/9*16;
pbyLength=1;

if isfield(opts.export, pbaspect)
    pbxLength=pbaspect(1);
    pbyLength=pbaspect(2);
end

subplot(3+numPeriods,1,1);
if helperCoDep==1
    j=1;
    for i=1:length(plotBack)-1
        if i==opts.coPoint % Plot black line in background for specific point
            k=mod(i,ylen)-TmeStmp(1)+1;
            thickness = (TmeStmp(2)-TmeStmp(1))/500;
            xdif=x(k+1)-x(k);
            xp = [x(k)-thickness*xdif x(k)+thickness*xdif x(k)+thickness*xdif x(k)-thickness*xdif];
            yss = ymid(j);
            if isfield(opts.export, "yaxis") && opts.export.yaxis(1,2) ~= inf
                diff = 5*maxyheight(j);
            else
                diff = maxyheight(j);
            end
            yp = [yss-diff yss-diff yss+diff+part_factor*diff*2 yss+diff+part_factor*diff*2];
            patch(xp, yp, [0 0 0] , 'EdgeColor','none', 'FaceAlpha', 1)
        end
        maxc = max(abs(plotBack((j-1)*ylen+1:j*ylen)));
        if maxc > 0
            yss = ymid(j);
            diff = maxyheight(j);
            if plotBack(i)>0
                negCont = yss;
                posCont = yss+plotBack(i)/maxc*diff;
            else
                negCont = yss+plotBack(i)/maxc*diff;
                posCont = yss;
            end
            colInd = ceil(plotBack(i)/maxc*(floor(opts.discrNmb/2)))+(floor(opts.discrNmb/2)+1);
        else
            yss = ymid(j);
            negCont = yss;
            posCont = yss;
            colInd = (floor(opts.discrNmb/2)+1);
        end
        
        if mod(i,ylen)~=0 && mod(i,ylen)>=TmeStmp(1)
            k=mod(i,ylen)-TmeStmp(1)+1;
            xdif=x(k+1)-x(k);
            xp = [x(k)-1/2*xdif x(k)+1/2*xdif x(k)+1/2*xdif x(k)-1/2*xdif];
            yp = [negCont negCont posCont posCont];
            patch(xp, yp, colors(colInd,:) , 'EdgeColor',colors(colInd,:), 'FaceAlpha', 1)
        elseif mod(i,ylen)==0
            j=j+1;
            subplot(3+numPeriods,1,j)
        else
        end
    end
end

subplot(3+numPeriods,1,1)
hold on
opts.export.exportName = strcat(exportName, "_y");
if isfield(opts.export, "yaxis")
    y.lims(3:4) = opts.export.yaxis(1,:);
end
if isfield(opts.export, "dashedLines")
    plot(ones(length(y.mu),1)*opts.export.dashedLines(1,1), "LineWidth",opts.lineWidth*0.5, 'Color', [0.84,0.84,0.84], 'LineStyle','--')
    plot(ones(length(y.mu),1)*opts.export.dashedLines(1,2), "LineWidth",opts.lineWidth*0.5, 'Color', [0.84,0.84,0.84], 'LineStyle','--')
    tempExport = opts.export;
    tempExport.dashedLines = opts.export.dashedLines(1,:);
else
    tempExport = opts.export;
end
if length(opts.delta)>1 % PLOT given delta distribution on top of existing plot
    if isfield(opts.export, "yaxis")
        maxdel = (opts.export.yaxis(1,2)-opts.export.yaxis(1,1))*0.1;
        basevalue = opts.export.yaxis(1,2);
    else
        maxdel = (y.lims(4)-y.lims(3))*0.1;
        basevalue = y.lims(4);
    end
    y.lims(4) = basevalue+maxdel;
    plotthis = maxdel*(opts.delta-1)/max(opts.delta-1)+basevalue;
    plotthis = plotthis(TmeStmp(1):TmeStmp(2));
    area(plotthis,basevalue, 'FaceColor',[0.6 0.6 0.6], 'FaceAlpha',0.3, 'EdgeColor','none')
end
plot_dist(y, type,samplesColored = opts.samplesColored, nmbsamples=opts.nmbsamples,  ylim=y.lims(3:4), pbaspect=[pbxLength pbyLength 1], lineWidth=opts.lineWidth*1, export=tempExport)

if isfield(opts.export, "format") && opts.export.format == "FIG"; title("input data"); end

subplot(3+numPeriods,1,2)
hold on
opts.export.exportName = strcat(exportName, "_LT");
if isfield(opts.export, "yaxis")
    LT.lims(3:4) = opts.export.yaxis(2,:);
end

if isfield(opts.export, "dashedLines")
    plot(ones(length(y.mu),1)*opts.export.dashedLines(2,1), "LineWidth",opts.lineWidth*0.5, 'Color', [0.84,0.84,0.84], 'LineStyle','--')
    plot(ones(length(y.mu),1)*opts.export.dashedLines(2,2), "LineWidth",opts.lineWidth*0.5, 'Color', [0.84,0.84,0.84], 'LineStyle','--')
    tempExport = opts.export;
    tempExport.dashedLines = opts.export.dashedLines(2,:);
    plot_dist(LT, type, samplesColored = opts.samplesColored, nmbsamples=opts.nmbsamples, ylim=LT.lims(3:4), pbaspect=[pbxLength pbyLength 1], lineWidth=opts.lineWidth*1, export=tempExport)
else
    plot_dist(LT, type, samplesColored = opts.samplesColored, nmbsamples=opts.nmbsamples, ylim=LT.lims(3:4), pbaspect=[pbxLength pbyLength 1], lineWidth=opts.lineWidth*1, export=opts.export)
end
if isfield(opts.export, "format") && opts.export.format == "FIG"; title("trend component"); end

for i=1:numPeriods
    subplot(3+numPeriods,1,2+i)
    %nexttile(2+i)
    hold on
    temp.mu = ST.mu(:,i);
    temp.sigma = ST.sigma(:,i);
    temp.Sigma = ST.Sigma(:,:,i);
    if (type == "comb") || (type == "spaghetti")
    temp.samples = ST.samples(:,:,i);
    end
    opts.export.exportName = strcat(exportName, "_ST", string(i));
    if isfield(opts.export, "yaxis")
        ST.lims(3:4,i) = opts.export.yaxis(2+i,:);
    end
    if isfield(opts.export, "dashedLines")
        plot(ones(length(y.mu),1)*opts.export.dashedLines(2+i,1), "LineWidth",opts.lineWidth*0.5, 'Color', [0.84,0.84,0.84], 'LineStyle','--')
        plot(ones(length(y.mu),1)*opts.export.dashedLines(2+i,2), "LineWidth",opts.lineWidth*0.5, 'Color', [0.84,0.84,0.84], 'LineStyle','--')
        tempExport = opts.export;
        tempExport.dashedLines = opts.export.dashedLines(2+i,:);
        if ~isempty(opts.delta) && isfield(yLTSTR, "AHat")
            temp.sigmaNew = ST.sigmaNew(:,i);
        end
        plot_dist(temp, type,samplesColored = opts.samplesColored, nmbsamples=opts.nmbsamples, ylim=ST.lims(3:4,i),pbaspect=[pbxLength pbyLength 1], lineWidth=opts.lineWidth*1, export=tempExport)
    else
        if ~isempty(opts.delta) && isfield(yLTSTR, "AHat")
            temp.sigmaNew = ST.sigmaNew(:,i);
        end
        plot_dist(temp, type,samplesColored = opts.samplesColored, nmbsamples=opts.nmbsamples, ylim=ST.lims(3:4,i),pbaspect=[pbxLength pbyLength 1], lineWidth=opts.lineWidth*1, export=opts.export)
    end
    if isfield(opts.export, "format") && opts.export.format == "FIG"; title(strcat("seasonal component"," ", string(i))); end
end

subplot(3+numPeriods,1,3+numPeriods)
hold on
opts.export.exportName = strcat(exportName, "_R");
if isfield(opts.export, "yaxis")
    R.lims(3:4) = opts.export.yaxis(3+numPeriods,:);
end
if isfield(opts.export, "dashedLines")
    plot(ones(length(y.mu),1)*opts.export.dashedLines(3+numPeriods,1), "LineWidth",opts.lineWidth*0.5, 'Color', [0.84,0.84,0.84], 'LineStyle','--')
    plot(ones(length(y.mu),1)*opts.export.dashedLines(3+numPeriods,2), "LineWidth",opts.lineWidth*0.5, 'Color', [0.84,0.84,0.84], 'LineStyle','--')
    tempExport = opts.export;
    tempExport.dashedLines = opts.export.dashedLines(3+numPeriods,:);
    plot_dist(R, type, samplesColored = opts.samplesColored, nmbsamples=opts.nmbsamples, ylim=R.lims(3:4), pbaspect=[pbxLength pbyLength 1], lineWidth=opts.lineWidth*1, export=tempExport)
else
    plot_dist(R, type, samplesColored = opts.samplesColored, nmbsamples=opts.nmbsamples, ylim=R.lims(3:4), pbaspect=[pbxLength pbyLength 1], lineWidth=opts.lineWidth*1, export=opts.export)
end
if isfield(opts.export, "format") && opts.export.format == "FIG"; title("residual component"); end

if opts.export.export && isfield(opts.export, "format") && opts.export.format == "FIG"
    curFig = gcf;
    curFig.Position = [34 58 1566 915];
    exportgraphics(curFig, strcat(opts.export.exportPath , exportName,".png"),...
        'BackgroundColor','white', 'Resolution',300);
end

% Change exportName to original input opts.export.exportName
opts.export.exportName = exportName;
if opts.plotCov
    plotCov(yLTSTR, numPeriods, opts)
end

if opts.plotCor
    yLTSTR = plotCor(yLTSTR, numPeriods, opts);
end

if opts.plotCorLength
    plotCorLength(yLTSTR, numPeriods, opts)
    if opts.export.export && isfield(opts.export, "format") && opts.export.format == "FIG"
        curFig = gcf;
        curFig.Position = [34 58 1566 915];
        exportgraphics(curFig, strcat(opts.export.exportPath , opts.export.exportName,"_CorLength",".png"),...
            'BackgroundColor','white', 'Resolution',900);
    end
end
fprintf('##### \n')
fprintf('##### Plotting DONE\n')
fprintf('##### \n')
fprintf('################################################################################# \n')

end % Function


%%%%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotCov(yLTSTR, numPeriods, opts)
% This function PLOTCOV creates a plot of the covariance matrix
% yLTSTR.Sigma
nmbColors = opts.discrNmb; % max 128
colorLut = ColorLUT(1);
colorLut = colorLut(round(linspace(1,257,nmbColors)),:); %specifies color discretization
figure
if ~opts.export.export; sgtitle('Covariance Matrix'); end
imagesc(yLTSTR.Sigma)

% White lines between block matrices y, LT, STs, R
sublength = length(yLTSTR.mu)/(numPeriods+3);
hold on
for i=1:numPeriods + 2
    plot(((i-1)*sublength+0.5+sublength)*ones(length(yLTSTR.mu),1), ':', 'Color','white', 'LineWidth', opts.lineWidth)
    plot(((i-1)*sublength+sublength+0.5)*ones(length(yLTSTR.mu),1),linspace(1,length(yLTSTR.mu), length(yLTSTR.mu)), ':', 'Color','white', 'LineWidth', opts.lineWidth)
end

colormap(colorLut); 
max_el = max(max(abs(yLTSTR.Sigma)));
if abs(max_el)>0
    clim([-max_el max_el])
end
pbaspect([1 1 1])
colorbar

if opts.export.export
    if isfield(opts.export, "exportName")
        exportName = opts.export.exportName;
    else
        exportName = "plot_dist";
    end
    opts.export.exportName = strcat(exportName, "_CovPlot");
    opts.export.pbaspect = [1 1 1];
    exportPlot(opts)
end
end % Function plotCov

function yLTSTR = plotCor(yLTSTR, numPeriods, opts)
% This function PLOTCOR creates a plot of the correlation matrix using the
% covariance matrix yLTSTR.Sigma and the function corComp that transforms
% the covariance matrix to the correlation matrix.
% This is part of the correlation exploration for UASTL.
yLTSTR = corComp(yLTSTR);

nmbColors = opts.discrNmb; % max 128
colorLut = ColorLUT(1);
colorLut = colorLut(round(linspace(1,257,nmbColors)),:); %128

figure
if ~opts.export.export; sgtitle('Correlation Matrix'); end
imagesc(yLTSTR.corMat);

% White lines between block matrices y, LT, STs, R
sublength = length(yLTSTR.mu)/(numPeriods+3);
hold on
for i=1:numPeriods + 2
    plot(((i-1)*sublength+0.5+sublength)*ones(length(yLTSTR.mu),1), ':', 'Color','white', 'LineWidth', opts.lineWidth)
    plot(((i-1)*sublength+sublength+0.5)*ones(length(yLTSTR.mu),1),linspace(1,length(yLTSTR.mu), length(yLTSTR.mu)), ':', 'Color','white', 'LineWidth', opts.lineWidth)
end

colormap(colorLut); 
clim([-1 1])
pbaspect([1 1 1])
colorbar
if opts.export.export
    if isfield(opts.export, "exportName")
        exportName = opts.export.exportName;
    else
        exportName = "plot_dist";
    end
    opts.export.exportName = strcat(exportName, "_CorPlot");
    opts.export.pbaspect = [1 1 1];
    exportPlot(opts)
end
end % Function plotCor

function yLTSTR = corComp(yLTSTR)
% This function CORCOMP computes the correlation matrix given the
% covariance matrix yLTSTR.Sigma.
cInvMult = zeros(size(yLTSTR,1),1);
for k = 1:size(yLTSTR.Sigma,1)
    c = yLTSTR.Sigma(k,k); 
    if(c < 1e-12)
        cInvMult(k) = 0;
    else
        cInvMult(k) = 1/sqrt(c);
    end
end
yLTSTR.corMat = diag(cInvMult)*yLTSTR.Sigma*diag(cInvMult);
end % Function corComp

function plotCorLength(yLTSTR, numPeriods, opts)
% This function PLOTCORLENGTH creates a plot of the correlation length for
% each component of the correlation matrix.
% This is part of the correlation exploration for UASTL.

sbplotNmb = numPeriods + 3;
lngthData = floor(length(yLTSTR.mu)/sbplotNmb);
if ~isfield(yLTSTR,"corMat")
    yLTSTR = corComp(yLTSTR);
end

signCorMat = sign(yLTSTR.corMat);
maxdim = length(signCorMat);
FullCorLength = zeros(maxdim,1);
VertR = zeros(maxdim, 1);
VertL = zeros(maxdim, 1);
for i=1:maxdim
    curSignEl = signCorMat(i,i);
    if abs(curSignEl)>0
        vert1 = 1;
        vert2 = 1;
        decB = mod(i,lngthData);

        % We only look in the current Data/Trend/Seasonal/Residuum
        % block. The following filters the running variables.
        if decB == 0
            varR = 0;
            varL = lngthData;
        elseif (decB < lngthData/2) || (decB > lngthData/2)
            varR = lngthData - decB;
            varL = decB;
        else
            varR = floor(lngthData/2);
            varL = floor(lngthData/2);
        end

        for j=1:varR-1 % look right
            if curSignEl == signCorMat(i,i+j)
                vert1 = vert1+1;
            else
                break;
            end
        end
        for j=1:varL-1 % look left
            if curSignEl == signCorMat(i,i-j)
                vert2 = vert2+1;
            else
                break;
            end
        end
        VertR(i) = vert1;
        VertL(i) = vert2;
        FullCorLength(i) = vert1+vert2-1;
    end
end

nmbColors = floor(opts.discrNmb/2)+1; % max 128
out = cell(length(yLTSTR.corMat),1);
for i=1:length(yLTSTR.corMat)
    ColInd = NaN(lngthData,1);
    k=2; 
    if (VertL(i)>=1) || (VertR(i)>=1) 
        for j = VertL(i)-1:-1:1
            ColInd(k) = min(ceil(yLTSTR.corMat(i,i-j)*nmbColors),nmbColors);
            k=k+1;
        end
        ColInd(k) = min(ceil(yLTSTR.corMat(i,i)*nmbColors),nmbColors);
        k=k+1;

        for j=1:VertR(i)-1
            ColInd(k) = min(ceil(yLTSTR.corMat(i,i+j)*nmbColors),nmbColors);
            k=k+1;
        end
        ColInd = ColInd(~isnan(ColInd));
        out{i} = ColInd;
    end
end
plotmat = zeros(length(yLTSTR.corMat),lngthData);
ysize=zeros(length(out),1);
for i=1:length(out)
    if ~isempty(out{i})
        insert = out{i};
        ysize(i)=length(insert);
        plotmat(i,1:length(insert)) = insert;
    end
end

plotmat = flipud(plotmat');

colors = ColorLUT(1);
colors = colors(129:257,:); % only one-sided ColorLUT
colors = colors(round(linspace(1,128,nmbColors)),:);
colors = [1 1 1; colors];

figure
if ~opts.export.export; sgtitle('Correlation Length'); end
for i=1:sbplotNmb
    subplot(sbplotNmb, 1, i)
    maxI = max(ysize( opts.stamp_indices((i-1)*opts.StmpLength+1:i*opts.StmpLength) ));
    image2plot = plotmat( end-maxI+1:end , opts.stamp_indices((i-1)*opts.StmpLength+1:i*opts.StmpLength) );

    imagesc(image2plot)
    hold on
    colormap(colors)
    
    VertLs = [VertL( opts.stamp_indices((i-1)*opts.StmpLength+1:i*opts.StmpLength) ), ...
        VertL( opts.stamp_indices((i-1)*opts.StmpLength+1:i*opts.StmpLength) )];
    VertLs = reshape(VertLs', [2*size(VertLs,1),1]);
    xs = [linspace(1,opts.StmpLength+1,opts.StmpLength+1);linspace(1,opts.StmpLength+1,opts.StmpLength+1)];
    xs = reshape(xs, [2*(opts.StmpLength+1),1]);
    xs([1, end])=[];
    hold on
    plot(xs-0.5, maxI+1-VertLs, "LineWidth",2, 'Color', 'black', 'LineStyle','-')
    if isfield(opts.export, "format") && opts.export.format == "FIG"
        if i==1
            title("input data"); 
        elseif i==2
            title("trend component"); 
        elseif i==sbplotNmb
            title("residual component"); 
        else
            title(strcat("seasonal component"," ", string(i-2))); 
        end
    end
end
if opts.export.export
    if isfield(opts.export, "exportName")
        exportName = opts.export.exportName;
    else
        exportName = "plot_dist";
    end
    for i=1:sbplotNmb
        subplot(sbplotNmb,1,i)
        opts.export.exportName = strcat(exportName, "_CorLength_", string(i));
        exportSubPlot(opts)
    end
end
end % Function plotCorLength

function exportSubPlot(opts)
% This function EXPORTSUBPLOT uses the current axis handle and exports it
% using the specified export options (opts.export.*).
    if ~isfield(opts.export, "format"); opts.export.format = "SUBs"; end
    ax = gca;
    set(ax,'XTick',[], 'XColor', [1 1 1])
    set(ax,'YTick',[], 'YColor', [1 1 1])
    set(ax, 'Color', 'none')
    set(ax,'Visible','on')
    if isfield(opts.export, "pbaspect")
        pbaspect(opts.export.pbaspect)
    else
        pbaspect([1 1 1])
    end
    colorbar off
    if opts.export.format == "SUBs"
        if opts.export.export
            exportgraphics(ax, strcat(opts.export.exportPath,opts.export.exportName,".pdf"),...
                'ContentType','vector', 'BackgroundColor','none');
        end
    elseif opts.export.format == "FIG"
        set(ax,'XColor', [0 0 0])
        set(ax, 'YColor', [0 0 0])
        box on
    end
end

function exportPlot(opts)
% This function EXPORTPLOT uses the current figure handle and exports it
% using the specified export options (opts.export.*).
ax = gca;
set(ax,'XTick',[], 'XColor', [1 1 1])
set(ax,'YTick',[], 'YColor', [1 1 1])
set(ax, 'Color', 'none')
set(ax,'Visible','on')
if isfield(opts.export, "pbaspect")
    pbaspect(opts.export.pbaspect)
else
    pbaspect([1 1 1])
end
colorbar off
if opts.export.export
    exportgraphics(ax, strcat(opts.export.exportPath,opts.export.exportName,".pdf"),...
        'ContentType','vector', 'BackgroundColor','none');
end
end % Function exportPlot

function map = ColorLUT(number) 
% This function COLORLUT imports the used color scheme generated via the
% code by Andy Stein and based on Kenneth Moreland's code for creating 
% diverging colormaps. It uses the 'CoolWarmFloat257.csv' accessible via 
% https://www.kennethmoreland.com/color-maps/CoolWarmFloat257.csv
A=readmatrix('CoolWarmFloat257.csv');
if number == 1
    map = diverging_map(A(:,1),[0.230, 0.299, 0.754],[0.706, 0.016, 0.150]);    % blue-white-red
end
end %Function ColorLUT