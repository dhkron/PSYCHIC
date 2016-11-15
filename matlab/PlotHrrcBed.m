% Plots hierarchys within the range
% @s & @e are in pixel units
% filter == true -> print only tads (filters out merges if true)
function [] = PlotHrrcBed(fMatrix,fHrrcBed,s,e,res,filter,figPath)
	hasFig = ( exist('figPath','var') && numel(figPath)>0 );

	if ischar(fMatrix)
		a = load(fMatrix);
	else
		a = fMatrix;
	end
	if ischar(fHrrcBed)
		bed = importdata(fHrrcBed);
	else
		bed = fHrrcBed;
	end
	if ischar(s)
		s = str2num(s);
	end
	if ischar(e)
		e = str2num(e);
	end
	if ischar(res)
		res = str2num(res);
	end
	if ischar(filter)
		filter = str2num(filter);
	end

	a_quality = bed.data(:);
	a_t = bed.textdata(:,4);
	a_s = bed.textdata(:,2);
	a_e = bed.textdata(:,3);

	is_tad = zeros(1,numel(a_t));
	for i = 1:numel(a_t)
		t = cell2mat(a_t(i));
		if t(1)=='T'
			is_tad(i) = 1;
		end
	end

	offset = 1-s; %So that s is at 1
	
	DisplayHeatmap(log2(1+a),[],s:e,'red');

	for i = 1:numel(a_s)
		s1 = floor(str2num(cell2mat(a_s(i)))/res)+1;
		e1 = floor(str2num(cell2mat(a_e(i)))/res)+1;
		if is_tad(i) || ~filter
			%Assumes: s1<e1
			if (s1 > s) && (e1 < e)
				%Plot square
				xp = [s1+offset e1+offset e1+offset NaN];
	        	        yp = [s1+offset s1+offset e1+offset NaN];
        	        	patch(xp,yp,'magenta','FaceAlpha',0.0,'EdgeColor','black','LineWidth',1.5,'EdgeAlpha',1);
			elseif (s1 > s) && (s1 < e)
				xp = [s1+offset e+offset e+offset NaN];
	        	        yp = [s1+offset s1+offset e+offset NaN];
        	        	patch(xp,yp,'magenta','FaceAlpha',0.0,'EdgeColor','black','LineWidth',1.5,'EdgeAlpha',1);
			elseif (e1 < e) && (e1 > s)
				xp = [s+offset e1+offset e1+offset NaN];
	        	        yp = [s+offset s+offset e1+offset NaN];
        	        	patch(xp,yp,'magenta','FaceAlpha',0.0,'EdgeColor','black','LineWidth',1.5,'EdgeAlpha',1);
			end
		end
	end

	title(sprintf('range %d-%d',res*(s-1)+1,res*e));
	g = gcf;
	if hasFig
		SaveFigure(g,figPath)
		close(gcf);
	end
