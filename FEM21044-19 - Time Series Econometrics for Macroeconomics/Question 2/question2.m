%Clear the workspace

clearvars all

% Load data
load('data.mat')
y = table2array(data(1:186,2:5))';
y_demeaned = y - mean(y,2).*ones(4,length(y(1,:)));

%% Question 2a

% Gecontroleerd met Anne
% 
%     for i = 1:4
% 
%         startingvalues=[1;1/3*var(y_demeaned(i,:));2/3*var(y_demeaned(i,:))];
% 
%         clearvars options
%         options  =  optimset('fmincon');
%         options = optimset(options, 'MaxFunEvals', 1e+6);
%         options = optimset(options, 'MaxIter', 1e+6);
%         options = optimset(options, 'TolFun', 1e-6);
%         options = optimset(options, 'TolX', 1e-6);
%         lb = [-10;0;0];
%         ub = [10;10;10];
%         [ML_parameters,ML_LogL]=fmincon('NegativeLogLikelihoodA', startingvalues,[],[],[],[],lb,ub,[],options,y_demeaned(i,:))
% 
%         format short
%         ML_std = sqrt( diag ( inv(  fdhess6('NegativeLogLikelihoodA',ML_parameters,y_demeaned(i,:)) )));
% 
%         [xi,P,~] = KalmanFilterA(ML_parameters, y_demeaned(i,:));
% 
%         prediction_1(:,i) = xi(1:185)*ML_parameters(1);
% 
%         prediction_4(:,i) = xi(1:182)*(ML_parameters(1)^4);
% 
%         correlation_1(:,i) = corr(prediction_1(:,i),y_demeaned(i,2:186)');
%         correlation_4(:,i) = corr(prediction_4(:,i),y_demeaned(i,5:186)');
% 
%     end
% 
%     correlation_1
%     correlation_4

%% Question 2b

% Gecontroleerd met Anne
% 
% startingvalues = [0.8;0.5;0.5;0.5;0.5;1;1;2;2;8;1;1;1];
% 
% lb = [-100;-100;-100;-100;-100;-100;-100;-100;-100;0;0;0;0];
% ub = [100;100;100;100;100;100;100;100;100;100;100;100;100];
% 
% [params, logl] = optimize('NegativeLogLikelihoodB',startingvalues,y_demeaned,[lb,ub]);

%% Question 2c

% Gecontroleerd met Anne

% Select the parameters with the lowest (or highest maximized) likelihood

% [~,index] = min(logl);
% optimized_parameters = params(:,index);
% 
% [~,P,~] = KalmanFilterB(optimized_parameters,y_demeaned);
% P(:,:,186)

%% Question 2d

%% Question 2d

% %from 2b: Get params and logL(parameters estimated with ML, log Likelihood)
% %from 2c: Select the parameters with the lowest (or highest maximized) likelihood
% 
% %from fminuncon, fmincon and fminsearch
% [~,index] = min(logl);
% optimized_parameters = params(:,index); 
% %get xi for optimized_parameters uising KF
% [xi,P,~] = KalmanFilterB(optimized_parameters,y_demeaned);
% 
% %set F and H using optimized parameters
% F=zeros(5,5);
% F(1,1)=optimized_parameters(1);
% F(2,2)=optimized_parameters(2);
% F(3,3)=optimized_parameters(3);
% F(4,4)=optimized_parameters(4);
% F(5,5)=optimized_parameters(5);
% 
% H=zeros(5,4);
% H(1,1)=optimized_parameters(6);
% H(1,2)=optimized_parameters(7);
% H(1,3)=optimized_parameters(8);
% H(1,4)=optimized_parameters(9);
% H(2,1)=1;
% H(3,2)=1;
% H(4,3)=1;
% H(5,4)=1;
% 
% %make one step ahead prediction for y_demeaned, 187
% xi_1step_b=F*xi(:,186);
% y_demeaned_187=transpose(H)*xi_1step_b;
% 
% %load the data to include 187th observation
% y_data = table2array(data(:,2:5))';
% %demean met 1:186 mean
% y_demeaned_data_subst = y_data - mean(y,2).*ones(4,length(y_data(1,:)));
% 
% %substitute 1 step ahead forecasted value of y4 in y_demeaned_e
% y_demeaned_data_subst(4,187)=y_demeaned_187(4);
% 
% %estimate parameters with new data 1:187
% 
% startingvalues = [0.8;0.5;0.5;0.5;0.5;1;1;2;2;8;1;1;1];
% 
% lb = [-100;-100;-100;-100;-100;-100;-100;-100;-100;0;0;0;0];
% ub = [100;100;100;100;100;100;100;100;100;100;100;100;100];
% 
% [params, logl] = optimize('NegativeLogLikelihoodB',startingvalues,y_demeaned_data_subst,[lb,ub]);
% [~,index] = min(logl);
% ML_parameters_d = params(:,index); 
% 
% [xi_d_subst,~]=KalmanFilterB(ML_parameters_d, y_demeaned_data_subst);
% 
% %voorspel y187 met nieuwe parameters
% F_d=zeros(5,5);
% F_d(1,1)=ML_parameters_d(1);
% F_d(2,2)=ML_parameters_d(2);
% F_d(3,3)=ML_parameters_d(3);
% F_d(4,4)=ML_parameters_d(4);
% F_d(5,5)=ML_parameters_d(5);
% 
% H_d=zeros(5,4);
% H_d(1,1)=ML_parameters_d(6);
% H_d(1,2)=ML_parameters_d(7);
% H_d(1,3)=ML_parameters_d(8);
% H_d(1,4)=ML_parameters_d(9);
% H_d(2,1)=1;
% H_d(3,2)=1;
% H_d(4,3)=1;
% H_d(5,4)=1;
% 
% xi_1step_subst_data=F_d*xi_d_subst(:,186);
% xi_1step_subst_data_2 = H'*xi_d_subst(:,187);
% y_demeaned_187_d=transpose(H_d)*xi_1step_subst_data;
% y_187_d= y_demeaned_187_d(4) + mean(y(4,:)); 
% 
% difference_d=y_demeaned_data_subst(:,187)-y_demeaned_187_d;

%% Question 2e

startingvalues = [0.8;0.5;0.5;0.5;0.5;1;1;2;2;8;1;1;1;1;1;1;1];

lb = [-100;-100;-100;-100;-100;-100;-100;-100;-100;0;0;0;0;0;0;0;0];
ub = [100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100;100];

[params, logl] = optimize('NegativeLogLikelihoodE',startingvalues,y_demeaned,[lb,ub]);

%Select the parameters with the lowest (or highest maximized) likelihood
[~,index] = min(logl);
optimized_parameters = params(:,index);

f = [optimized_parameters(1,1), optimized_parameters(2,1), optimized_parameters(3,1), optimized_parameters(4,1), optimized_parameters(5,1)];
h = [optimized_parameters(6,1), optimized_parameters(7,1), optimized_parameters(8,1), optimized_parameters(9,1)];
H = [h;eye(4)]; 
F = diag(f);

[xi,~] = KalmanFilterE(optimized_parameters, y_demeaned);

prediction_1 = H'*F*xi(:,1:185);
prediction_4 = H'*(F^4)*xi(:,1:182);    

for s = 1:4
    
   correlation_1(:,s) = corr(prediction_1(s,:)',y_demeaned(s,2:186)'); 
   correlation_4(:,s) = corr(prediction_4(s,:)',y_demeaned(s,5:186)');
   
end

%%

correlation_1
correlation_4

%% Question 2f
% 
% F = diag([0.85;0.85;0.85;0.85]) + ones(4,4)*1/10
% Q = 2*cov(y_demeaned')
% R = 4*cov(y_demeaned')
% 
% %F = 0.9*ones(4,4) + 0.1 * rand()*ones(4,4); % persistence of the state
% %Q = (1/2 + rand() ) * 1/3 * cov(y_demeaned'); % vaguely inspired by the data
% %R = (1/2 + rand() ) * 2/3 * cov(y_demeaned'); % vaguely inspired by the data
% 
% for k=1:1000       
%       [F,Q,R]=EM_step(F,Q,R,y_demeaned);
% end
% 
% [xi,~] = KalmanFilterF(F,Q,R,y_demeaned);
% 
% H = eye(4);
% 
% prediction_1 = H'*F*xi(:,1:185);
% prediction_4 = H'*(F^4)*xi(:,1:182);
% 
% for s = 1:4
%     
%    correlation_1(:,s) = corr(prediction_1(s,:)',y_demeaned(s,2:186)'); 
%    correlation_4(:,s) = corr(prediction_4(s,:)',y_demeaned(s,5:186)');
%    
% end
% 
% correlation_1'
% correlation_4'
% 
% F
% Q
% R
% 
% format long g
% 
% NegativeLogLikelihoodF(F,Q,R,y_demeaned)