function [U dU] = derivchol(A, dA)
% function [U dU] = derivchol(A, dA)
%
% INPUTS:
%   A is symmetric definite positive
%   dA(:,:,k) are symmetric for all k=1,...,size(dA)
% OUTPUTS:
%   U: cholesky of A
%   dU: derivative in the direction dA
%
% Cholesky decomposition and its derivative. Beside comuting the Cholestky
% matrix U (i.e. A = U'*U, U upper triangular) this function compute also
% dU the derivative of the Cholesky decomposition of the square matrix A
% in the given direction dA

% % Select method
% if nargin<3 || isempty(method)
%     method = 0;
% end
% 
% if method==1
%     
%     % SOLVING LINEAR EQUATION SLOW
%     
%     %% Preparation step
%     % we can spend time here, it's is irrelevant because it does not need the
%     % value of input matrix A and dA
%     n = size(A,1);
%     % Transposed-matrix permutation
%     P = bsxfun(@plus,1:n,(0:n-1)'*n);
%     P = P(:);
%     % Indexes arrays needed for sparse generation later
%     [I J] = ndgrid(1:n,1:n);
%     % concatenation of n arrays: [1] [1 2] ... [1...n]
%     iup = nonzeros(triu(I));
%     % concatenation of n arrays: [0] [1 1] ... [n-1...n-1]
%     l = nonzeros(triu(J));
%     n2 = n*n;
%     m = (n+n2)/2;
%     l = reshape(l-1,[1 m]);
%     I = bsxfun(@plus,I(:,iup),l*n);
%     J = repmat(1:m,n,1);
%     subs = logical(triu(ones(n)));
%     
%     %% Direct method
%     % Solving U'*X + X'*U = dA
%     %   U is trianguar superior cholesky decompition of A
%     %   X is triagular superior
%     
%     % This part is speed critical
%     tic
%     U = chol(A);
%     U = U.';
%     u = U(:,iup);
%     % LL is matrix such that LL*vector(triu(X)) = U'*X
%     LL = sparse(I(:),J(:),u(:),n2,m);
%     % LR is matrix such that LR*vector(triu(X)) = X'*U
%     LR = sparse(P(I(:)),J(:),u(:),n2,m);
%     % large inversion here, matrix size is (n*n) x m, which is about
%     % n^2 x n^2
%     x = (LL+LR) \ dA(:);
%     % Assign the triangular part
%     dU = zeros(n);
%     dU(subs) = x;
%     
% else
    % derivative by differientation Cholesky algorithm

% Robust Cholesky's factorization    
s = 0;
while true
    [U r] = chol(A);
    % preventive treatment, usually we do not need it
    if r==0
        break
    end
    try
        sm = eigs(A,1,'sa');
    catch %#ok
        sm = sqrt(eps(norm(A)));
    end
    se = eps(max([diag(A); norm(A)]));
    s = max([s sm se]);
    A = A + s*speye(size(A));
    s = s*2;
end

if nargout>=2 % derivative is requested
    cls = class(A);
    n = size(A,1);
    p = size(dA,3);
    dU = zeros([n n p],cls);
    % Loop on the matrices of derivative
    for k=1:p
        dUk = zeros([n n],cls);
        dAk = dA(:,:,k);
        % We duplicate the Cholesky's loop here
        for i=1:n
            subs = 1:i-1;
            j = i:n;
            ds = dAk(i,j)-(dUk(subs,i).'*U(subs,j)+U(subs,i).'*dUk(subs,j));
            Uii = U(i,i);
            dUii = ds(1)/(2*Uii);
            dUk(i,i) = dUii;
            dUk(i,j(2:end)) = (ds(2:end) - dUii*U(i,j(2:end))) / Uii;
        end % i for-loop
        dU(:,:,k) = dUk;
    end % k for-loop
end % if nargout>=2

%end

end % derivchol