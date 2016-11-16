function N = qrnull(R)
% function N = qrnull(R)
% Find the null space of the row-full-rank matrix R by QR method
% R must have less rows then columns (underdetermined system)
% Author: Bruno Luong <brunoluong@yahoo.com>
% History: 
%          04-Dec-2009: original

[m n] = size(R);
try
    [q r p] = qr(R.'); %#ok
catch %#ok
    % Older matlab version (anterior to 2009B) do not support
    % sparse-QR with permutation
    [q r p] = qr(full(R).'); %#ok
end
N = q(:,end+(m-n+1:0));
