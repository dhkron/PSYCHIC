function varargout = DerivBKnotDeriv(knots, k, r, dknots, alpha, lambda)
% [Dr td kd dDr subs] = DerivBKnotDeriv(knots, k, r, dknots)
% [alphad td kd dalphad subs] = DerivBKnotDeriv(knots, k, r, dknots, alpha)
% [alphad td kd dalphad subs dlambda] = ...
%             DerivBKnotDeriv(knots, k, r, dknots, alpha, lambda)
% Take the r^th derivative of a spline function (with respect to the
% abscissa coordinates x) of braket t_{j} and order k
%
% INPUTS:
%   knots: vector, knots points, must be ascending sorted
%   k-1: "order" of the spline (k is scalar)
%   r: order of the derivative
%   dknots: directions of the knot derivative, ranged in p columns
%   alpha: vector of length n := length(knots)-k, optional coefficients of
%         the spline basis
%   lambda: (optional) dual variable vector of size n-r
% OUTPUTS:
%   Dr: sparse matrix (n-r) x n, derivative matrix that maps coefficients
%       the original coefficients to the derivative coefficients
%       Note that Dr class will be double regardless the type of the knots
%   td: subsknots that can be used later to evaluate the derivative
%   kd: order of the derivative
%   dDr: cell of knot-derivative, each kth element corresponds to a
%        derivative in a direction dknots(:,k)
%   subs: index of subsknots, i.e. td = t(subs)
%   dalphad: matrix (n-r) x p, derivative when right-multiply by alpha
%   dlambda: matrix (n x p), derivative when left-multiply by lambda
%
% If spline coefficients alpha is provides, the coefficients of the
% derivative alphad of the spline is returned, otherwise Dr is returned

%%
nknots = length(knots);
n = nknots-k;

%%
if r>=k
    % zeros sparse matrix
    %Dr = sparse([],[],[],n-r,n);
    error('DERIVB: order derivation r(=%d) must be < k(=%d)', r, k);
end

if nargin<4
    dknots = zeros([nknots 0],class(knots));
end

coefin = nargin>=5 && ~isempty(alpha);

% Reshape in rows
knots = reshape(knots,1,[]);

if isequal(size(dknots),[1 nknots])
    p = 1;
else
    p = size(dknots,2);
    if size(dknots,1)~=nknots
        error('DerivBKnotDeriv: dknots must have %d rows', nknots);
    end
    dknots = dknots.';
end
    
%% 
if p==0
    % Compute the matrix only, no knot-derivative
    Dr = speye(n);
    for nu=1:r
        ij = nu+1:n;
        dt = knots(ij+(k-nu))-knots(ij);
        h = (k-nu)./dt;
        d = n-nu;
        HL = spdiags(h(:)*[-1 1],[0 1],d,d+1);
        Dr = HL*Dr;
    end
    
    DDR = {};
    
else
    
    DDR = cell(1,p);
    In = speye(n,n);
    Zn = sparse(n,n);
    %
    % Loop on vector of directions
    %
    for v=1:p
        Dr = In;
        dDr = Zn;
        % Loop on order of the derivatives
        for nu=1:r
            
            ij = nu+1:n;
            dk = (knots(ij+(k-nu))-knots(ij)); % vector
            ddk = (dknots(v,ij+(k-nu))-dknots(v,ij)); % vector
            h = (k-nu)./dk; % vector
            dh = -h.*ddk./dk; % vector
            
            d = n-nu;
            i = [1:d; ...
                 1:d];
            j = [1:d; ...
                 2:d+1];
            
            HL = sparse(i, j, [-h; h], d, d+1); % matrix
            dHL = sparse(i, j, [-dh; dh], d, d+1); % matrix
            
            dDr = dHL*Dr + HL*dDr;
            Dr = HL*Dr;
        end  % for-loop on nu
        
        % Because the matrices are sparse, we need to store
        % each of them in a cell
        DDR{v} = dDr;
    end % for-loop on v
    
end

%%
subs = r+1:length(knots)-r;
td = knots(subs);
kd = k-r;

%%
if coefin
    alpha = alpha(:);
    alphad = Dr*alpha;
    
    cls = class(knots);
     
    dalphad = zeros([n-r p],cls);
    for v=1:p
        dalphad(:,v) = DDR{v}*alpha;
    end
    
    % Product with the dual variable lambda
    if nargin>=6 % lambda is provided
        lambda = reshape(lambda, n-r, []).'; % row
        nl = size(lambda,1);
        dlambda = zeros([n nl p],cls);
        for v=1:p
            dlambda(:,:,v) = (lambda*DDR{v}).';
        end        
        out = {alphad td kd dalphad subs dlambda};
    else
        out = {alphad td kd dalphad subs};
    end
    
    varargout = out;

else
    % Call with unique derivative direction, we return the matrix
    % (no cell is needed)
    if p==1
        DDR = DDR{1};
    end
    varargout = {Dr td kd DDR subs};
end

