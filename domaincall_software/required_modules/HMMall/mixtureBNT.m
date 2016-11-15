 load 'mixtureBNT.mat'
 size(walkingX)
 size(runningX)
 dag = [ 0 1 1 ; 0 0 1 ; 0 0 0 ];
 discrete_nodes = [1 2];
 nodes = [1 : 3];
 node_sizes=[ 2 2 31];
 bnet = mk_bnet(dag, node_sizes, 'discrete', discrete_nodes);
 bnet.CPD{1} = tabular_CPD(bnet,1);
 bnet.CPD{2} = tabular_CPD(bnet,2);
 bnet.CPD{3} = gaussian_CPD(bnet, 3);
 %bnet.CPD{3} = gaussian_CPD(bnet, 3,'cov_type','diag');
 trainingX = walkingX(1:100,:);
 trainingX(101:200,:)=runningX(1:100,:);
 trainingC(1:100) = 1;   %% Class 1 is walking
 trainingC(101:200) = 2; %% Class 2 is running
 testX(1:20,:) = walkingX(101:120,:);   %% The first 20 are walking
 testX(21:40,:) = runningX(101:120,:);  %% The next 20 are running
 training= cell(3,length(trainingX));
 training(3,:) = num2cell(trainingX',1);
 training(1,:) = num2cell(trainingC,1);
 engine = jtree_inf_engine(bnet);
 maxiter=10;     %% The number of iterations of EM (max)
 epsilon=1e-100; %% A very small stopping criterion
 [bnet2, ll, engine2] = learn_params_em(engine,training,maxiter,epsilon);
 class0= cell(3,1); %% Create an empty cell array for observations
 class1 = class0;
 class2 = class0;
 class1{1} = 1;     %% The class node is observed to be walking
 class2{1} = 2;     %% The class node is observed to be running
 for i=1:100
   sample1=sample_bnet(bnet2,'evidence',class1);
   sample2=sample_bnet(bnet2,'evidence',class2);
   modelX(i,:)=sample1{3}';
   modelX(i+100,:)=sample2{3}';
 end
 figure
 subplot(2,1,1);
 plot(trainingX);
 subplot(2,1,2);
 plot(modelX);
 evidence=class0;   %% Start out with nothing observed
 for i=1:40
   evidence{3}=testX(i,:)';
   [engine3, ll] = enter_evidence(engine2,evidence);
   marg = marginal_nodes(engine3,1);
   p(i,:)=marg.T';
 end
 figure;
 subplot(2,1,1);
 plot(testX);
 hold
 plot(p(:,1));  %% Plot the output of the walking classifier
 subplot(2,1,2);
 plot(testX);
 hold
 plot(p(:,2));  %% Plot the output of the running classifier
