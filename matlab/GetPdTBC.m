function [pr_tad, pr_bg pr_crs] = GetPdTBC(a_t,a_b,a_c,MIN_DIAG,MAX_DIAG)
pr_tad = zeros(1,numel(MIN_DIAG:MAX_DIAG));
pr_bg = zeros(1,numel(MIN_DIAG:MAX_DIAG));
pr_crs = zeros(1,numel(MIN_DIAG:MAX_DIAG));
for i = MIN_DIAG:MAX_DIAG
	%cleanNan, cleanZero. NO LOG
	%dgn = GetDiag(a,i,1,1,0,0);
	%norm = sum(dgn);
	dgn = GetDiag(a_t,i,1,1,0,0);
	pr_tad(i) = sum(dgn);
	dgn = GetDiag(a_b,i,1,1,0,0);
	pr_bg(i) = sum(dgn);
	dgn = GetDiag(a_c,i,1,1,0,0);
	pr_crs(i) = sum(dgn);

	norm = pr_tad(i) + pr_bg(i) + pr_crs(i);
	pr_tad(i) = pr_tad(i)/norm;
	pr_bg(i) = pr_bg(i)/norm;
	pr_crs(i) = pr_crs(i)/norm;
end

end
