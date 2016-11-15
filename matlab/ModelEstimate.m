function [] = ModelEstimate(matPath,bedPath,bedPathOut,emMatOut,res)

	if ischar(res)
		res = str2num(res);
	end

	[RES_Hrr,~,~] = RebuildMatrix(matPath,bedPath,0,bedPathOut,res);

	Log('Writing model estiamted matrix');
	dlmwrite(emMatOut,RES_Hrr);
	Log();

end
