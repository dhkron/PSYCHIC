function [] = ModelEstimate(matPath,bedPath,bedPathOut,emMatOut,res,bilinear)

	if ischar(res)
		res = str2num(res);
	end
	[RES_Hrr,~,~] = RebuildMatrix(matPath,bedPath,0,bedPathOut,res,bilinear);

	Log('Writing model estiamted matrix');
	dlmwrite(emMatOut,RES_Hrr);
	Log();

end
