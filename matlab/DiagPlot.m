% Returns an array of means along all diagonals
function [L,K] = DiagPlot(m)
	if size(m,1) == 1
		L = m;
		K = m;
		return
	elseif size(m,2) == 1
		L = m';
		K = m';
		return
	end
	N = size(m,2)+size(m,1)-1;
	L = NaN*ones(1,N);
	K = NaN*ones(1,N);
	D1 = size(m,1);
	D2 = size(m,2);
	for i=1:N
		dg = i-1 - (D1-1);
		L(i) = nanmean(diag(m,dg));
		K(i) = nanmedian(diag(m,dg));
	end
end
