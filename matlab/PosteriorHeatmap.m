function [a_lpr,a_pdt,a_pdb] = PosteriorHeatmap(a,MIN_DIAG,MAX_DIAG,probBg,meansBg,sigmaBg,probIn,meansIn,sigmaIn,noNorm)

a_size = size(a,1);
a_pdt = zeros(a_size);
a_pdb = zeros(a_size);
for i = MIN_DIAG:MAX_DIAG
	start_pos = a_size*i + 1;
	diag_elem = GetDiag(a,i,0,0,1,1); %Still maps zero to zero

	pdt_diag = normpdf(diag_elem,meansIn(i),sigmaIn(i))*probIn(i);
	pdt_diag(diag_elem==0) = 0;
	a_pdt(start_pos:a_size+1:end) = pdt_diag;

	pdb_diag = normpdf(diag_elem,meansBg(i),sigmaBg(i))*probBg(i); %without prob it is liklihood
	pdb_diag(diag_elem==0) = 0;
	a_pdb(start_pos:a_size+1:end) = pdb_diag;
end

%This has a major effect on the results, a positive one
if ~exist('noNorm','var') || noNorm == false
	a_pdsum = a_pdt + a_pdb; 
	a_pdt = a_pdt ./ a_pdsum;
	a_pdb = a_pdb ./ a_pdsum;
end
a_pdt(isnan(a_pdt)) = 0;
a_pdb(isnan(a_pdt)) = 0;
a_lpr = log2(a_pdt)-log2(a_pdb); %Why log2?
a_lpr(isnan(a_lpr)) = 0;

if exist('box','var')
	a_pdt = a_pdt(box,box);
	a_pdb = a_pdb(box,box);
	a_lpr = a_lpr(box,box);
end
end
