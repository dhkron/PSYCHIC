function Log(msg)
	if exist('msg','var')
		fprintf(msg);
		fprintf('... ');
		tic;
	else
		elapsed = toc;
		minutes = floor(elapsed/60);
		seconds = mod(elapsed,60);
		if minutes == 0
			elstr = sprintf('%2.4f seconds',seconds);
		elseif minutes == 1
			elstr = sprintf('1 minute %2.2f seconds',seconds);
		else
			elstr = sprintf('%d minutes %2.2f seconds',minutes,seconds);
		end
		fprintf('Done! %s\r\n',elstr)
	end
end
