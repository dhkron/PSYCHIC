function a_gmm = AnalyzeGMM(a,MIN_DIAG,MAX_DIAG,k)
a_gmm = {};
warning('off','stats:gmdistribution:FailedToConverge');
for i = MIN_DIAG:MAX_DIAG
	%dgn = diag(a,i);
	%dgn_clean = dgn(dgn>0); %Zeros are bad for fitting??
	%dgn_clean = log(dgn_clean); %Log-normal % MUST BE CONSISTENT HERE
	dgn_clean = GetDiag(a,i,1,1,1,1);%Rm Zeros&Nan, Log % BE CONSISTENT IN LOG NORMAL +1 or +0

	if numel(dgn_clean)>1

		options = statset('MaxIter',1000);
		try
			GMModel = fitgmdist(dgn_clean,k,'Options',options);
			a_gmm{i}=GMModel;
		catch
			a_gmm{i}.mu = mean(dgn_clean);
			a_gmm{i}.Sigma = var(dgn_clean);
		end

		[~,msgid] = lastwarn;
		if strcmp(msgid, 'stats:gmdistribution:FailedToConverge')
			fprintf('Failed to converge at diag %d\r\n',i);
			warning('')
		end
	else
		a_gmm{i}.mu = 0;
		a_gmm{i}.Sigma = 0;
	end
end
w = warning('query','last');
end
