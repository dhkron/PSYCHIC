function [B dB dalpha] = BernKnotDeriv(x, knots, j, k, dknots, alpha, lambda, ...
                                       leftflag) %#codegen
% [B dB] = BernKnotDeriv(x, knots, j, k, dknots)
% [s ds] = BernKnotDeriv(x, knots, j, k, dknots, alpha)
% [s ds dalpha] = BernKnotDeriv(x, knots, j, k, dknots, alpha, lambda)
%
% Compute Bernstein polynomial basis by de De Casteljau's algorithm
% and its derivative with respect to the knots
%
% INPUTS:
%   x: vectors/array, point coordinates at which the function is to be
%      evaluated
%   knots: vector, knots points, must be ascending sorted
%   j: vector, vector of spatial index, must be in [1:length(knots)-k]
%         if it's empty all the basis functions are computed.
%   k-1: "order" of the spline (k is scalar)
%      k==1 -> piecewise constant
%      k==2 -> linear
%      k==4 -> cubic
%   dknots: increment of knots, each column is the direction where the
%          derivative is computed. In order to compute the Jacobian
%          user must provide a basis vectors for DKNOT 
%   alpha: vectors of size n:=length(t)-k, optional coefficients of the
%          spline basis
%   lambda: dual variable, same number of elements as x
%
% OUTPUTS:
%   B: (m x n) where m is length(x), n is length(j)              
%       Each column of B is the basis function B_j,k
%   dB: the derivative of B with respect to the direction given by
%       dknots
%   s: spline function, same dimension as x
%   ds: the derivative of s with respcte to knots
%   dalpha: left product the derivative dB (of B) with the dual lambda,
%           dimension (n x p)
%
%   If t is not monotonically non-decreasing, output will be an empty arrays
%
%   Note: B_j,k has support on [knots(j),knots(j+k)[
%
%  Call with TRUE for 8th argument: B = BernKnotDeriv(..., TRUE) to compute
%  the "Left B-spline", which has support on ]knots(j),knots(j+k)].
%  The left B-spline is needed e.g. to evaluate recursion integral between
%  two B-splines.      
%
% Author: Bruno Luong
%   31-May-2010: correct bug for LEFT Bspline is called with scalar

% Check if the coefficients are provided by user
coefin = nargin>=6 && ~isempty(alpha);

if coefin
    if size(alpha,1)==1
        alpha = alpha.';
    end
end
    
% Class of x
cls = class(x);

%%
if isvector(x) %&& ~coefin
    szx = length(x);
else
    szx = size(x);
end
x = x(:);
m = size(x,1);

%%
% Max possible value of indice
maxj = length(knots)-k;
if isempty(j) || coefin
    % all the index j
    j = 1:maxj;
    js = j;
else
    js  = sort(j(:));
end

% Reshape in column vector
if isequal(size(dknots),[1 length(knots)])
    dknots = dknots(:);
end
if isempty(dknots)
    dknots = zeros(length(knots),0,cls);
end
p = size(dknots,2);
dknots = reshape(dknots,[1 size(dknots)]);

% left and right bracket
jmin = js(1);
jmax = js(end);
% Check
if jmin<1 || jmax>maxj
    error('BERNSTEIN: j must be within [%d,%d]', 1, maxj);
end

%%
% Spcial case, we ignore the Dirac distribution
if k<=0    
    if coefin     
        B = zeros([szx size(alpha,2)],cls);
        dB = zeros([szx size(alpha,2) p],cls);
    else
        B = zeros([szx numel(j)],cls);
        dB = zeros([szx numel(j) p],cls);
    end
    return
end

%%
% k=1: step functions (piecwise constant)

B = zeros(m,jmax+k-jmin,cls);

if nargin>=8 && leftflag
    % Left B-spline
    tt = knots(jmax+k:-1:jmin);
    if issorted(-tt)
        [trash col] = histc(-x,-tt);
        col = length(tt)-col; % Correct BUG, 31/05/2010
    else
        B = [];
        dB = [];
        return
    end
else
    % Default, right B-spline
    tt = knots(jmin:jmax+k);
    if issorted(tt)
        [trash col] = histc(x,tt); %#ok
    else
        B = [];
        dB = [];
        return
    end
end

row = find(col>=1 & col<=size(B,2));
col = col(row);
%B(sub2ind(size(B),row,col)) = 1;
B(row+(col-1)*m) = 1;

% Derivative vanishes for k==0
dB = zeros([size(B) p],cls);

% dB has three dimensions: 
% - first -> abscissa (x)
% - second -> basis (j)
% - third -> knot derivative

%%
%bsxops(1);
for kk=2:k
    % recursion
    for jj=jmin:jmax+k-kk
        dt = knots(jj+kk-1)-knots(jj); % scalar  
        if dt~=0
            dkleft = dknots(1,jj,:);
            ddt = dknots(1,jj+kk-1,:)-dkleft; % (1 x 1 x p)
            w1 = (x-knots(jj)) / dt; % (m x 1 x 1)
            % dw1 = (-dknots(1,jj,:) - ddt.*w1) / dt;  % (m x 1 x p)
            dw1 = bsxfun(@times, ddt, w1);
            dw1 = bsxfun(@minus, -dkleft, dw1) / dt;
        else
            w1 = zeros(size(x),cls);
            dw1 = zeros([size(w1) p],cls);
        end
        dt = knots(jj+kk)-knots(jj+1);        
        if dt~=0
            dkright = dknots(1,jj+kk,:);
            ddt = dkright-dknots(1,jj+1,:); % (1 x 1 x p)
            w2 = (knots(jj+kk)-x) / dt; % (m x 1 x 1)
            %dw2 = (dknots(1,jj+kk,:) - ddt.*w2) / dt; % (m x 1 x p)
            dw2 = bsxfun(@times, ddt, w2);
            dw2 = bsxfun(@minus, dkright, dw2) / dt;
        else
            w2 = zeros(size(x),cls);
            dw2 = zeros([size(w2) p],cls);
        end
        ij = jj-jmin+1;
%         dB(:,ij,:) = dw1.*B(:,ij) + w1.*dB(:,ij,:) + ...
%                      dw2.*B(:,ij+1) + w2.*dB(:,ij+1,:);
        Bij = B(:,ij);
        Bijp1 = B(:,ij+1);
        dB(:,ij,:) = (bsxfun(@times, dw1, Bij) + ...
                      bsxfun(@times, dw2, Bijp1)) + ...
                     (bsxfun(@times, w1, dB(:,ij,:)) + ...
                      bsxfun(@times, w2, dB(:,ij+1,:))) ;
        B(:,ij) = w1.*Bij + w2.*Bijp1;
    end
end
%bsxops(0);

%%
% Map to original vector j
%[tf loc] = ismemberc(j, jmin:jmax); %#ok
if length(j) ~= size(B,2)
    loc = ismembc2(j, jmin:jmax);
    B = B(:,loc);
    dB = dB(:,loc,:);
end

% Basis
n = numel(j);
B = reshape(B, [szx n]);
dB = reshape(dB, [szx n p]);

% Multiply with coefficients to get the spline function
if coefin
    alpha = alpha(:);
    
    B = reshape(B,[],n);
    B = B*alpha;
    B = reshape(B, [szx size(alpha,2)]);

    if nargin>=7 % && nargout>=3 % lambda is provided
        % left product with the dual
        lambda = lambda(:).'; % row-vector
        dalpha = lambda*reshape(dB, szx, []);
        dalpha = reshape(dalpha, [n p]);
    end
    
    % permute the coefficient dimension to the end
    d = length(szx)+1;
    ip = 1:max(ndims(dB),d);
    ip([d end]) = ip([end d]);
    dB = reshape(permute(dB, ip), [], length(alpha));
    dB = dB*alpha;
    dB = reshape(dB, [szx p]);    

end

end
