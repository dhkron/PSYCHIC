function varargout = DerivB(t, k, r, alpha, columns)
% [Dr td kd subs] = DerivB(t, k, r)
% [alphad td kd subs] = DerivB(t, k, r, alpha)
%
% Take the r^th derivative of a spline function (with respect to the
% abscissa coordinates x) of braket t_{j} and order k
%
% INPUTS:
%   t: vector, knots points, must be ascending sorted
%   k-1: "order" of the spline (k is scalar)
%   r: order of the derivative
%   alpha: vectors of size n=length(t)-k, optional coefficients of the basis
%   columns: return only specified columns of the derivative
%     (corresponding to free-knots indexes)
% OUTPUTS:
%   Dr: sparse matrix (n-r) x n, derivative matrix that maps coefficients
%       the original coefficients to the derivative coefficients
%       Dr has at most (r+1) non-zero elements per row.
%   td: subsknots that can be used later to evaluate the derivative
%   kd: order of the derivative
%   subs: index of subsknots, i.e. td = t(subs)
%
% If spline coefficients alpha is provides, the coefficients of the
% derivative alphad of the spline is returned, otherwise Dr is returned

n = length(t)-k;

%% Special case
if r>=k
    % zeros sparse matrix
    Dr = sparse([],[],[],n-r,n);
    %error('DERIVB: order derivation r(=%d) must be < k(=%d)', r, k);
else
    
    %%
    Dr = speye(n);
    for nu=1:r
        ij = nu+1:n;
        dt = t(ij+(k-nu))-t(ij);
        h = (k-nu)./dt;
        d = n-nu;
        HL = spdiags(h(:)*[-1 1],[0 1],d,d+1);
        Dr = HL*Dr;
    end
end

%%
subs = r+1:length(t)-r;
td = t(subs);
kd = k-r;

%%
if nargin>=4 && ~isempty(alpha)
    alphad = Dr*alpha;
    varargout = {alphad td kd subs};
elseif nargin>=5 && ~isempty(columns)
    Dr = Dr(:,columns);
    varargout = {Dr td kd subs};
else
    varargout = {Dr td kd subs};
end

