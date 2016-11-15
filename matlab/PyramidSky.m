function [supermap,pyrIsTad,skyAreBlue] = PyramidSky(a_t,a_b,h,m)
% L(i,j) = sum( llr(pyramid is tad) + llr(sky are background) )

clr = triu(tril(ones(size(a_t)),h),m);
a_tad = fliplr(a_t .* clr);
a_bg = fliplr(a_b .* clr);

%figure; imagesc(fliplr(a_tad)); axis equal; colorbar;
%figure; imagesc(fliplr(a_bg)); axis equal; colorbar;

%pyrm = LlrPyramid(a_tad,h);
%clmn = LlrColumns(a_tad,h);

if sum(size(a_b)) == 2
	skyAreBlue = 0;
else
	skyAreBlue = LlrColumns(a_bg,2*h)-LlrPyramid(a_bg,2*h);
end
if sum(size(a_t)) == 2
	pyrIsTad = 0;
else
	pyrIsTad = LlrPyramid(a_tad,2*h);
end

supermap = pyrIsTad + skyAreBlue;
%supermap = skyAreBlue;
%supermap = pyrIsTad;

% /Begin normalization
% The dynamic programming algoritm makes this redundent I think	
if 0
	a_size = size(a_t,1);
	for i=1:a_size
		for j=1:a_size
			j2 = a_size+1-j;
			%nrm = 0.5*(i+j2-a_size-1)^2+1;
			%nrm2 = h*abs(i+j2-a_size-1);
			%pyrIsTad(i,j) = pyrIsTad(i,j)/nrm;
			%skyAreBlue(i,j) = skyAreBlue(i,j)/nrm2;
			%supermap(i,j) = supermap(i,j)/(abs(i-j)+1);
		end
	end

	%supermap = pyrIsTad + skyAreBlue;
	%supermap(isnan(supermap)) = 0;
	%supermap(isinf(supermap)) = 0;
end
% /End of normalization
return

end

function pyrm = LlrPyramid(a,h)
	a_size = size(a,1);
	pyrm = zeros(a_size);

	%Fill pyramid
	for x = a_size:-1:1 %Row
		for y = (a_size+1-x):-1:max(1,a_size+1-x-h) %Col
			pyrm(x,y) = LlrPyramidHelper(a,x,y);
		end
	end

	function psum = LlrPyramidHelper(a,i,j)
		if i+j == a_size+1
			psum = a(i,j);
		else
			%Using the next column
			%Therefore, fill up from bottom to top
			%And from diag elements to i=1
			psum = pyrm(i,j+1)+sum(a(i:end,j));
		end
	end

	pyrm = fliplr(pyrm);
end

function clmn = LlrColumns(a,h)
a_size = size(a,1);
clmn = zeros(a_size);

%Columns
for x = a_size:-1:1 %Row
		for y = (a_size+1-x):-1:max(1,a_size+1-x-h) %Col
			clmn(x,y) = LlrColumnsHelper(a,x,y);
		end
	end
	function psum = LlrColumnsHelper(a,i,j)
	psum = 0;
	if i+j == a_size+1
			psum = sum(diag(a,a_size+1-2*i));
		else
			%Fill columns starting from bottom row
			psum = clmn(i,j+1)+sum(diag(a,2*j-a_size-1))+sum(diag(a,2*j-a_size));
		end
	end

	clmn = fliplr(clmn);
end

