function f = DisplayHeatmap(a,range,box,scheme)
	f = figure('units','normalized','outerposition',[0 0 1 1]);%,'visible','off');
	hasRange = ( exist('range','var') && numel(range)>1 );
	hasBox = ( exist('box','var') && numel(box)>1);
	hasScheme = exist('scheme','var') && numel(scheme)>0;
	if hasRange && hasBox
		imagesc(a(box,box),range);
	elseif hasRange
		imagesc(a,range);
	elseif hasBox
		imagesc(a(box,box));
	else
		imagesc(a);
	end
	if hasScheme
		if strcmpi(scheme,'orange')
			orange = [ones(256,1),0.8*(256:-1:1)'/256+0.2,(256:-1:1)'/256];
			colormap(orange);
		else
			red = [ones(256,1),(256:-1:1)'/256,(256:-1:1)'/256];
			colormap(red);
		end
	end
	axis square;
	axis tight;
	%axis equal;
	colorbar;
end
