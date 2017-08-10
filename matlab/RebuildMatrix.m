%
% RES is actually the restored or rebuilt matrix.
% It is NOT the difference between that data and the rebuilt matrix.
%
function [RES_Hrr,RES_Tad,RES_Dxn] = RebuildMatrix(matPath,bedPath,dxnPath,bedOutPath,res,bilinear_fit_a)

	addpath(genpath('BSFKFolder')); %%%%%%%%%%Give some errors

    global bilinear_fit
    if nargin < 6
        bilinear_fit = true;
    else
        if isstr(bilinear_fit_a)
            bilinear_fit = str2num(bilinear_fit_a);
        else
            bilinear_fit = bilinear_fit_a;
        end
    end

	if dxnPath == 0
		skipDxn = true;
		RES_Tad = 0;
		RES_Dxn = 0;
	else
		skipDxn = false;
	end

	Log('Loading');

	nij = load(matPath);
	fBed = fopen(bedPath,'r');
	[BED,C] = textscan(fBed,'chr%s\t%f\t%f\t%s\t%f',Inf);

	firstMergeIndex = find(cellfun(@(x) strcmp(x,'Merge1'), BED{4}, 'UniformOutput', 1));

	chrNum = BED{1}{1};
	assert(numel(unique(BED{1}))==1); %Should not have more then one chromosome here

	myBounds = floor( [ cell2mat(BED(2)) cell2mat(BED(3)) ]/res )+1;
	if ~skipDxn
		myBoundsNoTad = myBounds(1:firstMergeIndex-1,:);
		dxnBounds = floor( importdata(dxnPath)/res+1 );
	end

	Log();

	%Tads and merges
	Log('Calculating Residual #1');
	[RES_Hrr,prm_hrr,X1] = GetResidual(myBounds,nij,res);
	Log();

	if ~skipDxn
		Log('Calculating Residual #2');
		[RES_Tad,prm_tad,X2] = GetResidual(myBoundsNoTad,nij,res);
		Log();

		Log('Calculating Residual #3');
		[RES_Dxn,prm_dxn,X3] = GetResidual(dxnBounds,nij,res);
		Log();
	end

	%RMSE
	fprintf('RMSE of Tads hierarchies %g\r\n',GetRMSE(RES_Hrr,nij,res));
	if ~skipDxn
		fprintf('RMSE of Tads without hierarchies %g\r\n',GetRMSE(RES_Tad,nij,res));
		fprintf('RMSE of Dixon with no hierarchies %g\r\n',GetRMSE(RES_Dxn,nij,res));
	end

	%Write beds
	BedWrite([bedOutPath],prm_hrr,chrNum,firstMergeIndex);
	if ~skipDxn
		BedWrite([bedOutPath '_2.bed'],prm_tad,chrNum,inf);
		BedWrite([bedOutPath '_3.bed'],prm_dxn,chrNum,inf);
	end
end

function [RMSE] = GetRMSE(RES,nij,res) %Changed by Dror Moran on 21.2.16
	clear_diag = CreateDiagMatrix(size(RES,1), 1, floor(5000000/res));
	RES(~clear_diag) = NaN; %Set anything above 5m to NaN
    % make sure there are no NaNs between 2nd diagonal and max distance
    assert(~any(isnan(RES(clear_diag == 1))));
    r1 = log2(nij(nij > 0));
    r2 = log2(RES(nij > 0));
	RES_sq = (r1 - r2).^2;
	RMSE = sqrt(nanmean(RES_sq(:)));
end

function [RES,Params,X] = GetResidual(bounds,nij,res)
	RES = NaN*ones(size(nij));
	Params = NaN*ones(9,size(bounds,1)+1); % +1 for sky
	for i=1:size(bounds, 1)
		X = bounds(i,1):bounds(i,2);
		if length(X)==1
			RES(X,X) = nij(X,X);
			tmp = [NaN;NaN;NaN;NaN;NaN;NaN;NaN];
		else
         	[newRES,tmp] = calcTadGradient(nij(X,X), RES(X,X), res);
            RES(X,X) = newRES;
		end
		Params(:,i) = [res*(bounds(i,1)-1); res*(bounds(i,2)-1); tmp];
	end
	%Sky
    fIdx = 1;
    lIdx = length(nij);
	X=fIdx:lIdx;
    [RES(X,X), tmp] = calcTadGradient(nij(X,X), RES(X,X), res);
	Params(:,end) = [res*(fIdx-1); res*(lIdx-1); tmp];
end

function [newRES, grad] = calcTadGradient(DAT, RES, res)
    global bilinear_fit
    % Don't take into account distances larger than 5Mb
    N = min(length(RES)-1, floor(5000000 / res));
	YY=ones(N,1);
	i = 1;
	count = 0;
	sumObj = zeros(1,2);
	L=NaN*ones(1,N);
	while i <= N
		% for each diagonal
		allDiag = diag(DAT,i);
		noRes = diag(RES, i);
		D=allDiag(find(isnan(noRes)));
		
		count = count + length(D);
		sumObj(1) = sumObj(1) + length(D)*i;
		sumObj(2) = sumObj(2) + nansum(D);
		if count >= 10 && sumObj(2) ~= 0
			L(round(sumObj(1)/count)) = sumObj(2)/count; 
			count = 0;
			sumObj = zeros(1,2);
		end
		i = i + 1;
    end
    if count > 0 && sumObj(2) ~= 0
		L(round(sumObj(1)/count)) = sumObj(2)/count; 
    end
    
	XX=res*[1:N]; 
	II=find(~isnan(L));
   
	%find non-NAN elements
	% regress
    
	pbreak =NaN;
	ValueCount = length(II);
	if ValueCount > 1
        if bilinear_fit
            %Split
            logXX = log2(XX(II));
            logL = log2(L(II));
            [pp , ~]=BSFK(logXX, logL, 2, 2);
            pbreak = pp.breaks(2);
		
            logXX1 = logXX(logXX < pbreak);
            logL1 = logL(logXX < pbreak);
            logXX2 = logXX(logXX >= pbreak);
            logL2 = logL(logXX >= pbreak);
            
       
            %Regression
            if length(logXX2) > 1 && length(logXX1) > 1 
                disp('bilinear');
                YY1 = ones(length(logL1),1);
                [b1,~,~,~,stat1]=regress(logL1', [logXX1', YY1]);
                YY2 = ones(length(logL2),1);
                [b2,~,~,~,stat2]=regress(logL2', [logXX2', YY2]);
                if b2(1) >= 0 || b1(1) >= 0
                    YY = ones(length(logL),1);
                    [b1,~,~,~,stat1]=regress(logL', [logXX', YY]);
                    b2=b1;
                    stat2=stat1;
                end
            else
                disp('linear');
                YY = ones(length(logL),1);
                [b1,~,~,~,stat1]=regress(logL', [logXX', YY]);
                b2=[NaN;NaN];
                stat2=NaN;
            end
        else
	        % linear fit
            [b1,stat1] = find_power_law(DAT, RES, res, N);
            b2 = [NaN;NaN];
            stat2 = NaN;
        end 
	elseif ValueCount == 1
		 b1 = [NaN;NaN]; b2 = [NaN;NaN]; stat1 = NaN; stat2 = NaN; display('TAD with less then 10 values!!!'); %%%%Maybe Change if happend
	else
		b1 = [NaN;NaN]; b2 = [NaN;NaN]; stat1 = NaN; stat2 = NaN;
	end
	
	%Rebuilt matrix
	newRES = zeros(length(RES));
	
	if ValueCount > 1
		b = [];
		for i=1:N
			oldDiag = diag(RES, i);
			newDiag = zeros(1, length(oldDiag));
			
			if log2(i*res) < pbreak || all(isnan(b2))
				b = b1;
			else
				b = b2;
			end
			newDiag(find(isnan(oldDiag))) = 2^(log2(i*res)*b(1)+b(2));
            newRES = newRES + diag(newDiag, -i) + diag(newDiag, i);
		end
    end
    assert(all(all(isnan(RES) | (newRES == 0) | (RES == newRES))));
	newRES = bsxfun(@max, RES, newRES);%????
	%The diagonal is not part of the model and therefore should be NaN, Dror 1.5.16 
	newRES(logical(eye(size(newRES)))) = NaN;
	grad = [b1; stat1(1); b2; stat2(1); pbreak];
end

function BedWrite(outPath,bed,chrNum,mergeIndex)
	%Add 'chr1' at col1, type at col4, than normal
	f = fopen(outPath,'w');
	for i = 1:size(bed,2)
		if i >= mergeIndex
			txt = 'Merge';
		else
			txt = 'TAD';
		end
		if i == size(bed,2)
			txt = 'Sky';
		end
		ln = bed(:,i);
		fprintf(f,'chr%s\t%d\t%d\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',chrNum,ln(1),ln(2),txt,ln(3),ln(4),ln(5),ln(6),ln(7),ln(8),ln(9));
	end
	fclose(f);
end

function [grad,stats] = find_power_law(DAT, RES, res, N)
    A = [];
    b = []; 
    for i=1:N
        for j=i+1:N
            if isnan(RES(i,j)) && DAT(i,j) > 0
                A = [A ; [log2(res * (j-i)) 1]];
                b = [b ; log2(DAT(i,j))];
            end
        end
    end
    if size(A, 1) == 0
        disp('EMPTY A!');
        stats = NaN;
        grad = [0 0]';
    else
        [grad,~,~,~,stats] = A\b;
    end
    
end
