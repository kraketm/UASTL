%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Uncertainty-Aware Seasonal-Trend Decomposition Based on Loess %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the main demo for 
% 'Uncertainty-Aware Seasonal-Trend Decomposition Based on Loess'.
%
% 
% This demo uses a simple, artificial (discrete) Gaussian process with 
% moderate uncertainty.

clear; close all; clc;
addpath('VIS/')


%% Input
rng('default')

% Model the mean mu
N = 200;
t = (1:N)';
trend = t/10;
periodic = 10*sin(2*pi*t/100);
noise = 2*(rand(N,1) - 0.5);
mu = trend + periodic + noise;

% Model the covariance matrix Sigma
sigma2 = 20*ones(N,1);
sigma = sqrt(sigma2);
Sigma = zeros(N,N);
ExQuKernel = @(x,y,sigma_i,sigma_j,l) sigma_i*sigma_j*exp(-0.5*norm(x-y)^2 / l^2);
for i = 1:N
    for j = 1:N
        Sigma(i,j) = ExQuKernel(t(i),t(j),sigma(i),sigma(j),5);
    end
end

% Create Gaussian distributed variable F
X.t = t;
X.mu = mu;
X.Sigma = Sigma;


%% UASTL
yLTSTR = uastl(X,100);


%% UASTL and Correlation exploration
export.export = false;                       % set to true if you want to export
export.exportPath = "figures/";
export.exportName = "fig_1";
export.pbaspect = [7 1 1];
export.format = "FIG";                      % set to "SUBs" for individual subplot export
% export.yaxis          - matrix for y-axis limits, e.g., export.dashedLines =[-20 20;0 40;-20 20;-7.5 7.5];
% export.dashedLines    - matrix to insert dashed lines for spefific y-values, e.g., export.dashedLines =[-15 15;5 15;-15 15;-5 5];

yLTSTR = plot_distributionmtx(yLTSTR, ...   % consisting of .mu and .Sigma
    1,'comb', ...                           % number of Periods and plot_type (comb/isoband/spaghetti)
    samplesColored=true, ...                % spaghettis colored logical
    nmbsamples=3, ...                       % nmb of spaghetti samples
    plotCov=false, ...                      % covariance matrix plot?
    plotCor=true, ...                      % correlation matrix plot?
    plotCorLength=false, ...                % correlation length plot?
    coPoint=100, ...                        % interactive point for correlation exploration
    export=export, ...                      % export set as above to export individual subplots
    lineWidth=2.5, ...                      % lineWidth Factor for all plots
    discrNmb=13 ...                         % binning of colorbar    
    );


%% UASTL and Sensitivity Analyis
% For the sensitivity analysis one needs to compute the second output of
% "uastl" once, which is the global matrix 'AHatGlobal'
[yLTSTR, AHatGlobal] = uastl(X, 100);
yLTSTR.AHat = AHatGlobal(:,1:length(yLTSTR.mu)/(1+3)); % to boost efficiency, only parts of the global matrix AHatGlobal are used

export.export = false;                       % set to true if you want to export
export.exportPath = "figures/";
export.exportName = "fig_2";
export.pbaspect = [7 1 1];
export.format = "FIG";                      % set to "SUBs" for individual subplot export
% export.yaxis          - matrix for y-axis limits, e.g., export.dashedLines =[-20 20;0 40;-20 20;-7.5 7.5];
% export.dashedLines    - matrix to insert dashed lines for spefific y-values, e.g., export.dashedLines =[-15 15;5 15;-15 15;-5 5];

% Create kernel for sensitivity analysis
pos = 100;
ampl = 1;
bdwidth = 20;
delta = ones(length(yLTSTR.mu)/(1+3),1);
delta(pos-bdwidth:pos+bdwidth) = 1 + ampl * exp(-0.5*abs( (1:2*bdwidth+1) - (bdwidth+1) ).^2 / (bdwidth/3)^2);
 
yLTSTR = plot_distributionmtx(yLTSTR, ...   % consisting of .mu and .Sigma
    1,'comb', ...                           % number of Periods and plot_type (comb/isoband/spaghetti)
    samplesColored=true, ...                % spaghettis colored logical
    nmbsamples=3, ...                       % nmb of spaghetti samples
    plotCov=false, ...                      % covariance matrix plot?
    plotCor=false, ...                      % correlation matrix plot?
    plotCorLength=false, ...                % correlation length plot?
    coPoint=pos, ...                        % interactive point for sensitivity analysis
    export=export, ...                      % export set as above to export individual subplots
    lineWidth=2.5, ...                      % lineWidth Factor for all plots
    discrNmb=13, ...                        % binning of colorbar
    delta=delta, ...                        % kernel for the sensitivity analysis
    timeStamp=[50,150] ...                  % shows only the specific interval
    );

