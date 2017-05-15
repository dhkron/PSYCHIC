function SaveFigure(handle, path)
	%handle.Units = 'centimeters';
	%handle.PaperUnits = 'centimeters';
	%handle.PaperPosition = handle.Position;

	[~,fname,fext] = fileparts(path);
	
	Log(['Saving file "',fname,fext,'"']);
	saveas(handle,path); %In matlab, the extension matters
	Log();

	%if findstr(path,'.png') && exist('export_fig')
	%	higherQ = regexprep(path,'\.png$','.better.png');
	%	Log('Saving file in better quality');
	%	export_fig(handle,'-a4',higherQ);
	%	Log();
	%end
end
