function [R dR dt] = DGD(t, k, dt, freeidx)
% [R dR dt] = DGD(t, k, freeidx) OR
% [R dR dt] = DGD(t, k, [], dt)
%
% Compute the simplified discrete Gram diagonal matrix (DGD) of B-spline
% basis of knots t
%
% INPUTS:
%   t: array of knots, must be ascending sorted (not checked)
%      elements could be repeated
%   k-1: the order of the spline
%   freeidx: index of (free) knots where the derivative of the Gram matrix
%            will be computed
%   dt: matrix of derivative directions (freeidx argument will be ignored)
%       size(dt,1) = size(t,1) (number of knots)
%       size(dt,2) = size(dR,3)
% OUTPUTS:
%   R: the DGD matrix of the basis, each diaginal term
%   dR: the matrix of derivative; arranged by columns of dt
%       size(dR,3) = size(dt,2)
%       size(dR(:,:,i)) = size(R) for all i

t = reshape(t,[],1);
nknots = size(t,1);
cls = class(t);

if k<1 || k>=nknots
    error('DGD: k(=%d) must be positive and smaller than %d', k, nknots)
end

if nargin >=4
    % Default direction from freeidx
    r = max(k-1,0);
    if isempty(freeidx)
        freeidx = 1+r:nknots-r;
    else
        if any(freeidx<1+r | freeidx>nknots-r)
            error('DGD: freeidx must comprise in (%d,%d)', 1+r, nknots-r);
        end
    end
    % vector of direction for derivative
    dt = zeros(nknots, numel(freeidx), cls);
    dt(freeidx + (0:numel(freeidx)-1)*nknots) = 1;
elseif nargin<3 || isempty(dt)
    dt = zeros(nknots,0);
elseif size(dt,1)~=nknots
    error('DGD: dt must have %d rows', nknots);
end

p = size(dt,2);

n = nknots-k;
lsupport = t(1+k:nknots) - t(1:nknots-k);
lsupport = lsupport / k;
R = zeros([n n], cls);
R(1:n+1:end) = lsupport;

dR = zeros([n*n p], cls);
dlsupport = dt(1+k:nknots,:) - dt(1:nknots-k,:);
dlsupport = dlsupport / k;
dR(1:n+1:end,:) = dlsupport;
dR = reshape(dR, [n n p]);

end % DGD

