% **********************************************
% ** WARNING: THIS CODE IS BASED ON TADTREE.m **
% **********************************************
%
% This file merges TADs based on the regression model
%
% It looks on the shared region between the TADs,
% and compare it against the background model
% The shared region should be something between a bgmodel and the two TADs around it
%
% Here I describe what I actually do.
% The magic function here is the mighty "RegressMergerHelper"
% 
% Input: TADs 
%	@bnd boundary locations
%	@a is the hic map
%	@box_s start of region of interest (inclusive)
%	@box_e end of region of interest (inclusive)
%	@res resolution of HiC map
%	
% Output: Hierarchy
%	@chrNum & @fBedPath are for output
%
% Examples: chr10 gettaEMT
% RegressMerger([11,38,54,85,100,120,139],a,151,300,40000,bgModel,9999999);
%
function [] = RegressMerger(bnd, a, box_s, box_e, res, bgModel, chrNum, fBedPath, doFig)
	hasFig = ( exist('doFig','var') && doFig );

	% Bound relative, always having 0 as start and 'end' as end
	% So giving relative bounds is possible, as long as start and end are well
	bounds = [0 bnd box_e-box_s+1];
	% The matrix @a is referenced by offset.
	% Boundaries are relative, so to access real data we add the offset
	offset = box_s-1;

	if hasFig	
		DisplayHeatmap(log2(1+a),0,box_s:box_e,'red');
		textDraw = {};
	end

	if exist('fBedPath','var') && numel(fBedPath)>0
		fBed = fopen(fBedPath,'w');
	else
		fBed = fopen('/dev/null','w'); % <3 Linux 
	end

	%fprintf('Range 1\t| Range2\t| V_DynProg\t| V_Super\n')
	%fprintf('-----------------------------------------------\n')
	import java.util.LinkedList;
	bl = LinkedList();

	%Generate the nice table of possible merges
	bound_matrix = zeros(numel(bounds)-2,5);
	bounds3 = [bounds(1:end-2) ; bounds(2:end-1) ; bounds(3:end) ];
	for i = 1:size(bounds3,2) 
		b = bounds3(:,i);
		s1 = b(1)+1;
		e1 = b(2);
		s2 = b(2)+1;
		e2 = b(3);
		%The merge parameter: Value of alpha of the merge area
		alpha_val = RegressMergeHelper(s1,e1,s2,e2,a,offset,bgModel);%
		bound_matrix(i,:) = [s1 ,e1, s2, e2, alpha_val];

		if hasFig	
			xp = [s1 e1 e1 NaN];
			yp = [s1 s1 e1 NaN];
			patch(xp,yp,'magenta','FaceAlpha',0.0,'EdgeColor','black','LineWidth',1.5,'EdgeAlpha',1);
			textDraw{end+1} = [e1 s1 0];
		end

		% offset-1 : 1 is mapped to 0, 2 is mapped to RES, etc.
		% With base of 501, 1 is actually 501, so add 500. Then remove 1 cause base 1 should be 0.
		% Print the bedfile
		M = sprintf('chr%s\t%d\t%d\tTAD\t%g\n',chrNum,(s1+offset-1)*res,(e1+offset-1)*res,0);
		fprintf(fBed,M);
	end
	%Add last TAD. First one is included in the loop. This is lazy coding.
	if i>0
		if hasFig	
			xp = [s2 e2 e2 NaN];
			yp = [s2 s2 e2 NaN];
			patch(xp,yp,'magenta','FaceAlpha',0.0,'EdgeColor','black','LineWidth',1.5,'EdgeAlpha',1);
			textDraw{end+1} = [e2 s2 0];
		end

		M = sprintf('chr%s\t%d\t%d\tTAD\t%g\n',chrNum,(s2+offset-1)*res,(e2+offset-1)*res,0);
		fprintf(fBed,M);
	end
	bound_matrix;

	mergeNumber = 1;
	while size(bound_matrix,1)>1
		sz = size(bound_matrix,1);
		max_alpha = 0;
		max_i = 0;
		max_alpha_s = -1; %Holds the start of the maximal row
		max_alpha_e = -1; %Holds the end of the maximal row
		%First, find the minimal alpha
		for i = 1:sz
			alpha_val = bound_matrix(i,5);
			
			if alpha_val > max_alpha
				max_alpha = alpha_val;
				max_i = i;
				max_alpha_s = bound_matrix(i,1);
				max_alpha_e = bound_matrix(i,4);
			end
		end
		if max_alpha == 0
			break
		end
		
		%Update table entries. Remember, each raw represents an INTERACTION
		% therefore, a merge is just the deletion of the interaction row
		%
		% 1. Update left, if has left
		% 2. Update right, if has right
		% 3. Remove row
		%
		%There are two possible steps, but sometimes only one is applicible
		% Do update left (upstream tad)
		if max_i > 1
			c = bound_matrix(max_i-1,:);
			
			s_p1 = max_alpha_s + offset;
			s_p2 = max_alpha_e + offset;
			s_p3 = c(1) + offset;

			new_alpha= RegressMergeHelper(c(1),c(2),c(3),max_alpha_e,a,offset,bgModel);
			bound_matrix(max_i-1,:) = [c(1),c(2),c(3),max_alpha_e,new_alpha];
		end
		% Do update right (downstream tad)
		if max_i < sz
			c = bound_matrix(max_i+1,:);

			s_p1 = max_alpha_s + offset;
			s_p2 = max_alpha_e + offset;
			s_p3 = c(4) + offset;
			
			new_alpha= RegressMergeHelper(max_alpha_s,c(2),c(3),c(4),a,offset,bgModel);%
			bound_matrix(max_i+1,:) = [max_alpha_s,c(2),c(3),c(4),new_alpha];
		end

		if hasFig
			xp = [max_alpha_s max_alpha_e max_alpha_e NaN];%max_alpha_e];
			yp = [max_alpha_s max_alpha_s max_alpha_e NaN];%max_alpha_s];
			patch(xp,yp,'magenta','FaceAlpha',0.0,'EdgeColor','black','LineWidth',1.5,'EdgeAlpha',1);
			textDraw{end+1} = [max_alpha_e, max_alpha_s, mergeNumber];	
		end

		%Debug
		fprintf('%d\t%d-%d\t%g\n',mergeNumber,max_alpha_s,max_alpha_e,max_alpha);
		
		bound_matrix = removerows(bound_matrix,'ind',max_i);
		
		M = sprintf('chr%s\t%d\t%d\tMerge%d\t%g\n',chrNum,(max_alpha_s+offset-1)*res,(max_alpha_e+offset-1)*res,mergeNumber,max_alpha);
		fprintf(fBed,M);

		mergeNumber = mergeNumber + 1;
	end

	if hasFig
		for c = textDraw
			drw = c{1};
			max_alpha_e = drw(1);
			max_alpha_s = drw(2);
			mergeNumber = drw(3);
			if mergeNumber == 0
				tsize = 1;
			elseif mergeNumber > 9
				tsize = 8;
			else
				tsize = 10;
			end
			t = text('position',[max_alpha_e, max_alpha_s],'fontsize',tsize,'string',sprintf('%d',mergeNumber),'HorizontalAlignment','center','Color','black','FontName','arial');
			tExt = t.Extent;
			rad = ceil(max([tExt(3), tExt(4)]))+1;
			rectangle('Position',[max_alpha_e-rad/2 max_alpha_s-rad/2 rad rad],'Curvature',[1,1,],'FaceColor','white');
			t = text('position',[max_alpha_e, max_alpha_s],'fontsize',tsize,'string',sprintf('%d',mergeNumber),'HorizontalAlignment','center','Color','black','FontName','arial');
		end
	end
end

function alpha_val = RegressMergeHelper(s1,e1,s2,e2,a,offset,bg)
	box1 = (s1+offset):(e1+offset);
	box2 = (s2+offset):(e2+offset);

	if numel(box1) == 1 || numel(box2) == 1
		alpha_val = 1;
		return
	end
	
	[~,La1] = DiagPlot(a(box1,box1));
	[~,La2] = DiagPlot(a(box2,box2));
	La1 = La1( (numel(La1)+3)/2 : end ); %fix because of triangle, see DiagPlotter.m
	La2 = La2( (numel(La2)+3)/2 : end ); %fix because of triangle, see DiagPlotter.m

	a_temp = NaN*ones(e2-s1+1);
	a_temp(1:(e1-s1+1),1:(e1-s1+1)) = a(box1,box1);
	a_temp((s2-s1+1):(e2-s1+1),(s2-s1+1):(e2-s1+1)) = a(box2,box2);
	[~,La] = DiagPlot(a_temp);

	La = La( (numel(La)+3)/2 : end ); %fix because of triangle, see DiagPlotter.m
	%DisplayHeatmap(log2(a_temp+1),0,[],'red');

	t_merge = a(box1,box2);
	[~,Lm] = DiagPlot(t_merge);

	%l = min([numel(Lm),numel(La),numel(bg)]);
	%LA = La(1:l)';
	%LM = Lm(1:l)';
	%BG = bg(1:l)';
	%alpha_val = regress( LM-BG, LA-BG );
	
	l1 = min([numel(Lm),numel(La1),numel(bg)]);
	LA1 = La(1:l1)';
	LM1 = Lm(1:l1)';
	BG1 = bg(1:l1)';
	alpha_val1 = regress( LM1-BG1, LA1-BG1 );
	
	l2 = min([numel(Lm),numel(La2),numel(bg)]);
	LA2 = La(1:l2)';
	LM2 = Lm(1:l2)';
	BG2 = bg(1:l2)';
	alpha_val2 = regress( LM2-BG2, LA2-BG2 );

	alpha_val = (alpha_val1*alpha_val2)^(0.5);
	if alpha_val > 1
		alpha_val = 2-alpha_val;
	end
end
