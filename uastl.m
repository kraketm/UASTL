function [XHat,AHatGlobal] = uastl(X,p,opts)
% This function UASTL computes uncertainty-aware seasonal-trend
% decomposition based on loess for Gaussian distributed data.
% The theory is described in the related manuscript "Uncertainty-Aware
% Seasonal-Trend Decomposition Based on Loess".
%
%
% [Xhat,AHatGlobal] = UASTL(X,p,opts)
%
% Input:
%   * X: uncertain data, Gaussian distributed variable with
%       - X.mu: mean vector
%       - X.Sigma: covariance matrix
%   * p: periods p_1,...,p_L
%   * options:
%       - robust: include robustness weights
%       - n_o: number of outer loops
%       - n_i: number of inner loops
%       - n_s: span for seasonal trends
%       - n_l: span for low-pass filter
%       - n_t: span for long-term trend
%       - postSmoothingSeasonal 
%       - postSmoothingSeasonal_n
%       - postSmoothingTrend
%       - postSmoothingTrend_n
%
% Output:
%   * XHat: 
%       - XHat.mu: mean vector
%       - XHat.Sigma: covariance matrix
%   * AHatGlobal
%


arguments
    X (1,1) struct
    p (1,:) double {mustBeInteger}
    opts.robust (1,1) logical = false
    opts.n_o    (1,1) double {mustBeInteger} = 1
    opts.n_i    (1,1) double {mustBeInteger} = 2
    opts.n_s    (1,:) double {mustBeInteger} = max(2*floor(ceil(length(X.mu)./p)/2)+1,5)
    opts.n_l    (1,:) double {mustBeInteger} = max(2*floor(p/2)+1,5)
    opts.n_t    (1,1) double {mustBeInteger} = 2*floor(((1.5*p(end))/(1 - 1.5/max(ceil(length(X.mu)/p(end)),5)))/2)+1;
    opts.postSmoothingSeasonal (1,1) logical = true
    opts.postSmoothingSeasonal_n (1,:) double {mustBeInteger} = max(2*floor((p/2)/2)+1,5)
    opts.postSmoothingTrend (1,1) logical = true
    opts.postSmoothingTrend_n (1,1) double {mustBeInteger} = max(2*floor((length(X.mu)/2)/2)+1,5)
end

fprintf('################################################################################# \n')
fprintf('##### Uncertainty-Aware Seasonal-Trend Decomposition Based on Loess (UASTL) ##### \n')
fprintf('################################################################################# \n')
fprintf('##### \n')
fprintf('##### Dimension:\t n = %i \n',length(X.mu))
fprintf('##### \n')
fprintf('##### Period %i:\t\t p_%i = %i \n',[1:numel(p);1:numel(p);p])
fprintf('##### \n')
fprintf('##### Robustness:\t %s \n',mat2str(opts.robust))
fprintf('##### \n')
fprintf('##### Outer Loop:\t n_o = %i \n',opts.n_o)
fprintf('##### Inner Loop:\t n_i = %i \n',opts.n_i)
fprintf('##### \n')
fprintf('################################################################################# \n')
fprintf('##### \n')
fprintf('##### UASTL START \n')
fprintf('##### \n')

if p ~= sort(p)
    error('UASTL: sort the periods from small to large')
end

if any([length(p) ~= length(opts.n_s), length(p) ~= length(opts.n_l), length(p) ~= length(opts.postSmoothingSeasonal_n)])
    error('UASTL: periods p, parameter n_s, parameter n_l, and parameter postSmoothingSeasonal_n must have the same size')
end

n = length(X.mu);
L = numel(p);

% Line 1: Initialize Embedding
nHat = (3+L)*n;
XHat.mu = zeros(nHat,1);
XHat.mu(1:n) = X.mu;
XHat.Sigma = zeros(nHat);
XHat.Sigma(1:n,1:n) = X.Sigma;

if nargout > 1
    AHatGlobal = eye(nHat);
end

% Line 2: Initialize weights
weights = ones(n,1);

% Line 3: Outer loop
for outer_loop = 1:opts.n_o
    fprintf('##### \t Outer Loop | iter: %i \n',outer_loop)

    % Line 4: Initialize STL Matrix AHat
    AHat = eye(nHat);

    % Line 5: Inner loop
    for inner_loop = 1:opts.n_i
        fprintf('##### \t\t Inner Loop | iter: %i ... ',inner_loop)
                
        % Line 6: Loop over periods
        for k = 1:L

            % Line 7: Update AHat regarding seasonal trends --- Steps 1-4
            
            % ------- Step 1: detrending()
            A_Deltat = zeros(3+L);
            A_Deltat(2+k,1) = 1;
            A_Deltat(2+k,2:2+k-1) = -1;
            A_Deltat(2+k,2+k) = 0;
        
            % ------- Step 2: cycle_subseries_smoothing(p_k,n_s,omega)
            E_ext = zeros(n+2*p(k),n); 
            E_ext(p(k)+1:p(k)+n,1:n) = eye(n);
            E_ext(1:p(k),1:p(k)) = eye(p(k));
            E_ext(end-p(k)+1:end,end-p(k)+1:end) = eye(p(k));

            E_split  = zeros(n+2*p(k));     
            B_ns = zeros(n+2*p(k));
            indx = 0;
            for i = 1:p(k)
                cycleSubseries = (i-p(k):p(k):i+floor((n-i)/p(k))*p(k)+p(k)) + p(k);
                lenCS = length(cycleSubseries);
                
                cycleSubseriesWeights = weights(i:p(k):end);
                B_ns(indx+1:indx+lenCS,indx+1:indx+lenCS) = loessmtx(lenCS,opts.n_s(k),2,cycleSubseriesWeights([1,1:end,end]));
                E_split(indx+1:indx+lenCS,cycleSubseries) = eye(lenCS);
                
                indx = indx+lenCS;
            end

            % ------- Step 3: cycle_subseries_low_pass_filtering(p_k,n_l)
            h = 3*(0:p(k))';
            h([1;end]) = h([1;end]) + [1;-2];
            h = [h; h(end-1:-1:1)] / (p(k)^2*3);
            A_L = convmtx(h,n);
            B_nl = loessmtx(n, opts.n_l(k),1);

            % -------  Step 4: cycle_subseries_detrending()
            P_1n = zeros(n,n+2*p(k));
            P_1n(:,p(k)+1:p(k)+n) = eye(n);

            A_p = (P_1n - B_nl*A_L')*E_split'*B_ns*E_split*E_ext;
        

            % update STL matrix AHut regarding seasonal trends
            A_S_id = eye(3+L);
            A_S_id(2+k,2+k) = 0;
            AHat = (kron(A_S_id,eye(n))+kron(A_Deltat,A_p)) * AHat;
 
        end % Loop over periods

        % Line 9: Update AHat regarding long-term trend --- Steps 5-6

        % ------- Step 5: deseasonalizing()
        A_tDash = zeros(3+L);
        A_tDash(2,1) = 1;
        A_tDash(2,3:3+L-1) = -1;
        A_tDash(2,2) = 0;

        % ------- Step 6: trend_smoothing(n_t,omega)
        A_T_id = eye(3+L);
        A_T_id(2,2) = 0;
        AHat = (kron(A_T_id,eye(n))+kron(A_tDash,loessmtx(n,opts.n_t,1,weights))) * AHat;

        fprintf('Done \n')
    end % Inner loop

    % Line 11: Update AHat regarding post smoothing of seasonal trends
    if opts.postSmoothingSeasonal
        for k = 1:L
            A_S_post = zeros(L+3);
            A_S_post(2+k,2+k) = 1;
            AHat = (kron(eye(L+3)-A_S_post,eye(n))+kron(A_S_post,loessmtx(n,opts.postSmoothingSeasonal_n(k),2))) * AHat;
        end
    end

    % Line 11: Update AHat regarding post smoothing of long-term trend
    if opts.postSmoothingTrend
        AHat = (kron(A_T_id,eye(n))+kron(A_tDash,loessmtx(n,opts.postSmoothingTrend_n,2))) * AHat;
    end

    % Line 12: Update AHat regarding residuum
    tmpmtx = eye(L+3);
    tmpmtx(L+3,1) = 1;
    tmpmtx(L+3,2:L+2) = -1;
    tmpmtx(L+3,L+3) = 0;
    AHat = kron(tmpmtx,eye(n)) * AHat; 

    % Line 13: Update embedding via linear transformation
    XHat.mu = AHat * XHat.mu;
    XHat.Sigma = AHat * XHat.Sigma * AHat';

    if nargout > 1
        AHatGlobal = AHat * AHatGlobal;
    end
   
    % Line 14: Update weights omega
    if opts.robust && outer_loop < opts.n_o
        r = mvnrnd(XHat.mu,XHat.Sigma,n);
        r = r(:,end-n+1:end);
        h = 6*median(abs(r),2);
        u = abs(r)./h;
        u2 = (1 - u.^2).^2;
        u2(u>1) = 0;
        weights = mean(u2)';
    end

end % Outer loop


fprintf('##### \n')
fprintf('##### UASTL DONE\n')
fprintf('##### \n')
fprintf('################################################################################# \n')


end % Function