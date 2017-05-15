% Diags from min to max
% To take the primary diagonal, set min=0
function clear_diag = CreateDiagMatrix(sizea, minD,maxD)
	clear_diag = triu(tril(ones(sizea),maxD),minD);
end
