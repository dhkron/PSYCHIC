function [f] = DebugDixonDomains(fMatrix,fDomains,res,box)
	if ischar(fMatrix) 
		a = load(fMatrix);
	else
		a = fMatrix;
	end
	if ischar(fDomains)
		b = load(fDomains);
	else
		b = fDomains;
	end
	hasBox = ( exist('box','var') && numel(box)>1);
	if ~hasBox
		box = 1:size(a,2);
	end
	f = DisplayHeatmap(log2(a+1),[0 15],box,'red');%Was orange for some reason
	hold on;
	ax = axis;
	if ischar(res)
		res = str2num(res);
	end
	for c = b'
		c = c/res;
		i = c(1);
		j = c(2);
		if i > box(1) && i < box(end)
			i = i-box(1)+1;
			plot([i,i],[ax(3),i],'g--');
			plot([i,ax(2)],[i,i],'g--');
		end
		if j > box(1) && j < box(end)
			j = j-box(1)+1;
			plot([j,j],[ax(3),j],'b-.');
			plot([j,ax(2)],[j,j],'b-.');
		end
	end
end
