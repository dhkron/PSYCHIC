function [dynmap, dynsol, dynmsk] = DynProgTAD(supermapraw,s,e,cost,c_dynmsk)
supermap = supermapraw(s:e,s:e);
dynmap = zeros(size(supermap));
dynsol = zeros(size(supermap));
dynmsk = zeros(1,size(supermap,1));%Will hold 1 for 'allowed' boundary solutions

supermapfl = fliplr(supermap);

if ~exist('c_dynmsk','var')
	c_dynmsk = ones(1,size(supermap,1));
end
%perms = {};
%for i=1:(e-s+1)
%	perms{i} = randperm(i);
%end

a_size = size(supermap,1);
for dg = 0:a_size
	for pos = 1:(a_size-dg)
		i = pos;%row
		j = dg+pos;%col

		if abs(i-j)<5
			dynmap(i,j) = supermap(i,j);
		else
			ptad = supermap(i,j);
			midtype = 0;
			%rng = i:j-1;
			%rng = rng(perms{numel(rng)});
			%for mid = rng
			for mid = i:j-1
				if c_dynmsk(mid) == 0
					continue;
				end
				flipj = a_size+1-j;
				v = dynmap(i,mid) + dynmap(mid+1,j) - cost;
				% + Interaction square value?
				% How does -100 affect the structure?
				%If I add a measure of the cross-tad interaction?
				%How can it find DISTANT interacting TADs?
				if v > ptad
					ptad = v;
					midtype = mid;
				end
			end
			dynmap(i,j) = ptad;
			dynsol(i,j) = midtype;
		end
	end
end
%figure; imagesc(dynmap); axis equal; colorbar;

disp('Finished! Printing solution...')
dynmsk = printSolution(dynsol,1,e-s+1,0); %offset=s-1
fprintf('\r\n');

end

function [dynmsk] = printSolution(dynsol,s,e,offset)
	dynmsk = zeros(1,e);%Will hold 1 for 'allowed' boundary solutions

	import java.util.LinkedList;
	q = LinkedList();
	q.add([s,e]);
	while q.size()>0
		item = q.remove();
		i = item(1);
		j = item(2);
		current = dynsol(i,j);
		if current == 0
			%Do printin and add to mask
			fprintf('%d-%d ',i+offset,j+offset)
			if abs(i-j) > 5
				if ~isempty(findall(0,'Type','Figure'))
					hold on;
					ax = axis;
					plot([i,i],[ax(3),i],'g--');
					plot([i,ax(2)],[i,i],'g--');
					plot([j,j],[ax(3),j],'b-.');
					plot([j,ax(2)],[j,j],'b-.');
					%xp = [i i j j];
					%yp = [i j j i];
					%patch(xp,yp,'magenta','FaceAlpha',0.1,'EdgeColor','none');
				end
			end
		else
			dynmsk(current) = 1;
			q.push([current+1 j]);
			q.push([i current]);
		end
	end
end
