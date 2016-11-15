function [a,a_llr,s] = FindNewTads(fMat,fSupersum,res,box,fNewDomains)
	hasBox = ( exist('box','var') && numel(box)>1 );

	Log('Loading files');
	a = load(fMat);
	s = load(fSupersum);
	Log();

	if ischar(res)
		res = str2num(res);
	end;

	if ~hasBox
		box = 1:size(a,2);
	end

	Log('Finding TADs');
	[d1 d2 d3]=DynProgTAD(s,box(1),box(end),0);
	Log();
	
	Log('Saving TADs');
	bounds = [1 find(d3) box(end)-box(1)+1];
	bounds2 = [bounds(1:end-1) ; bounds(2:end)];
	realBounds = (bounds2'-1)*res;
	dlmwrite(fNewDomains,realBounds,'precision',20);
	Log();
end

