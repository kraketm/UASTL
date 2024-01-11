function [B] = loessmtx(a,s,d,omega)
% This function LOESSMTX computes a matrix that performs the loess method.
% Loess (locally weighted scatterplot smoothing) smooths data based on 
% a local fitting with a linear or quadratic model.
% The theory is described in the related manuscript "Uncertainty-Aware
% Seasonal-Trend Decomposition Based on Loess".
%
%
% [B] = LOESSMTX(a,s,d,omega)
% 
% Input:
%   * a: related positions (if integer n is used: a = (1,2,...,n)^T)
%   * s: span
%   * d: degree of model (d=1 linear, d=2 quadratic)
%   * omega: external weights
%
% Output:
%   * B: loess matrix 
%


% 1. initialization and checks
if (isscalar(a) && floor(a)==a)
    n=a;
    a = (1:n)';
elseif isvector(a)
    a = a(:);
else
    error('loessmtx: related positions a must be a vector or an integer.')
end

n = length(a);

if nargin < 4
    omega = ones(n,1);
else
    omega = omega(:);
end

if s < 4
    error('loessmtx: span s is too small, it should be at least 4.')
end

s = floor(min(s,n));

B = zeros(n,n);

for i=1:n
    % 2. find the s-nearst neighbors of a_i
    [~,idx] = mink(abs(a-a(i)),s);
    
    % 3. center positions
    aTild = a(idx)-a(i);

    % 4. compute scaling
    aTildAbs = abs(aTild);
    scaling = (1 - (aTildAbs/max(aTildAbs)).^3).^3;

    % 5. define Vandermonde matrix
    V = aTild.^(0:d);

    % 6. compute diagonal matrix
    D = diag(scaling.*omega(idx));

    % 7. weighted linear least squares solution
    a_i = zeros(n,1);
    a_i(idx) = ((1:d+1)==1)*pinv(V'*D*V)*V'*D;

    % 8. insert row to loess matrix B
    B(i,:) = a_i;
end

end

