function dgn = GetDiag(a,i,cleanNan,clearZero,doLog,addOne)
	dgn = diag(a,i);
	if cleanNan
		dgn = dgn(~isnan(dgn));
	end	
	if clearZero
		dgn = dgn(dgn>0);
	end
	if doLog
		if addOne
			dgn = log(dgn+1); % This is ln a not log because: X lognormal => Y = lnX is normal
		else
			dgn = log(dgn);
		end
	end
end
