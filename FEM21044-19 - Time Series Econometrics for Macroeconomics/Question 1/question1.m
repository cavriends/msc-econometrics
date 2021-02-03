%% Load data

clearvars all

load('data.mat')
y = table2array(data(1:186,4)).';
dates = table2array(data(1:186,1)).';
dates2 = datenum('31/12/1899','dd/mm/yyyy')+dates(1:end)';

% Plotting
plot(dates2,y,'-o','MarkerSize',5,'MarkerEdgeColor','k')
dateFormat = 'mmmyy';
datetick('x',dateFormat);
startrow=datenum('1/1/1973','dd/mm/yyyy');
endrow=datenum('1/7/2019','dd/mm/yyyy');
axis([startrow,endrow,-inf,inf])
set(gca,'FontSize',10)
set(gca,'FontName','Times New Roman')

%% Question 1a

% First set of starting values
startingvalues = [0.8;0.8;0.95*mean(y);1.05*mean(y);std(y);std(y)];

lb = [0;0;-10;10;0;0];
ub = [1;1;20;20;20;20];

[ML_parameters_1,ML_LogL_1]=optimize('NegativeLogLikelihood',startingvalues,y,[lb,ub]);
% Select optimal parameters
[~,index] = min(ML_LogL_1);
optimal_1 = ML_parameters_1(:,index);

% Second set of starting values
startingvalues = [0.8;0.8;2*mean(y);-mean(y);1/2*std(y);2*std(y)];

[ML_parameters_2,ML_LogL_2]=optimize('NegativeLogLikelihood',startingvalues,y,[lb,ub]);
[~,index] = min(ML_LogL_2);
optimal_2 = ML_parameters_2(:,index);

%% Question 1b

p11   = optimal_1(1,1);
p22   = optimal_1(2,1);
mu    = [optimal_1(3,1);optimal_1(4,1)];
sigma = [optimal_1(5,1);optimal_1(6,1)];

[optimal_1(1,1);optimal_1(2,1);optimal_1(3,1);optimal_1(4,1);optimal_1(5,1);optimal_1(6,1)]

[smoothedxi_1] = Hamilton_smoother(p11,p22,mu,sigma,y);

% Plotting probability of being in state 2 for ML, startset 1
figure
plot(dates2,y,'-o','MarkerSize',5,'MarkerEdgeColor','k')
hold on
plot(dates2(1:end),smoothedxi_1(2,:),'r','LineWidth',2)
dateFormat = 'mmmyy';
datetick('x',dateFormat);
startrow=datenum('1/1/1973','dd/mm/yyyy');
endrow=datenum('1/7/2019','dd/mm/yyyy');
axis([startrow,endrow,-inf,inf])
set(gca,'FontSize',10)
set(gca,'FontName','Times New Roman')

p11   = optimal_2(1,1);
p22   = optimal_2(2,1);
mu    = [optimal_2(3,1);optimal_2(4,1)];
sigma = [optimal_2(5,1);optimal_2(6,1)];

[optimal_2(1,1);optimal_2(2,1);optimal_2(3,1);optimal_2(4,1);optimal_2(5,1);optimal_2(6,1)]

[smoothedxi_2] = Hamilton_smoother(p11,p22,mu,sigma,y);

% Plotting probability of being in state 2 for ML, startset 2
figure
plot(dates2,y,'-o','MarkerSize',5,'MarkerEdgeColor','k')
hold on
plot(dates2(1:end),smoothedxi_2(2,:),'r','LineWidth',2)
dateFormat = 'mmmyy';
datetick('x',dateFormat);
startrow=datenum('1/1/1973','dd/mm/yyyy');
endrow=datenum('1/7/2019','dd/mm/yyyy');
axis([startrow,endrow,-inf,inf])
set(gca,'FontSize',10)
set(gca,'FontName','Times New Roman')

% Plotting probability of being in state 2 for both startingpoints

figure
plot(dates2,y,'-o','MarkerSize',5,'MarkerEdgeColor','k')
hold on
plot(dates2(1:end),smoothedxi_1(1,:),'r','LineWidth',2)
hold on 
plot(dates2(1:end),smoothedxi_2(2,:),'g','LineWidth',2)
dateFormat = 'mmmyy';
datetick('x',dateFormat);
startrow=datenum('1/1/1973','dd/mm/yyyy');
endrow=datenum('1/7/2019','dd/mm/yyyy');
axis([startrow,endrow,-inf,inf])
set(gca,'FontSize',10)
set(gca,'FontName','Times New Roman')

%% Question 1c

% Set the starting values randomly

p11 = 0.8;
p22 = 0.8;
mu1 = 2*mean(y);
mu2 = -mean(y);
sigma1 = 1/2*std(y);
sigma2 = 2*std(y);

iteration = 0;

a= 1/2; %rand();
xi0_in=[a;1-a];

for k=1:10 % there is no k here, we just repeat this section 100 times
    % Increase counter
    iteration = iteration+1;
    % Run EM step
    [p11,p22,mu1,mu2,sigma1,sigma2,xi0_out]=EM_step(p11,p22,mu1,mu2,sigma1,sigma2,xi0_in,y);
    % Set in=out to make sure we update:
    xi0_in = xi0_out; 
    % Save the parameters
    parameters_EM(:,iteration) = [p11;p22;mu1;mu2;sigma1;sigma2;xi0_out];
    % print the latest parameters
    parameters_EM(:,end); 
    % Calculate the log-likelihood at these parameters
    LogL_EM(iteration) = - NegativeLogLikelihood2([p11;p22;mu1;mu2;sigma1;sigma2;xi0_out],y);
end

% Plot the last parameter estimate (EM)
parameters_EM(:,end)

% Check difference with ML (make sure use the same starting values)
optimal_2 - parameters_EM(1:6,end)

% Plot it
figure
plot(LogL_EM,'LineWidth',2) 

%figure
plot(parameters_EM(1:6,:)','LineWidth',2)

[filteredxi_1,~] = Hamilton_filter2(p11,p22,mu1,mu2,sigma1,sigma2,xi0_in,y);

% Plotting probability of being in state 2 for EM, startset 2
figure
plot(dates2,y,'-o','MarkerSize',5,'MarkerEdgeColor','k')
hold on
plot(dates2(1:end),filteredxi_1(2,:),'r','LineWidth',2)
dateFormat = 'mmmyy';
datetick('x',dateFormat);
startrow=datenum('1/1/1973','dd/mm/yyyy');
endrow=datenum('1/7/2019','dd/mm/yyyy');
axis([startrow,endrow,-inf,inf])
set(gca,'FontSize',10)
set(gca,'FontName','Times New Roman')

%% Question 1d

% Answer is in pdf

%% Question 1e

y = table2array(data(1:186,2:3))';

p11 = 0.8;
p22 = 0.8;
mu1 = 2*mean(y,2);
mu2 = -mean(y,2);
sigma1 = 1/2*cov(y');
sigma2 = 2*cov(y');

a= 1/2; %rand();
xi0_in=[a;1-a];

iteration = 0;

a= 1/2; %rand();
xi0_in=[a;1-a];
    
for k=1:1000 % there is no k here, we just repeat this section 100 times
    % Increase counter
    iteration = iteration+1;
    % Run EM step
    [p11,p22,mu1,mu2,sigma1,sigma2,xi0_out]=EM_step_E(p11,p22,mu1,mu2,sigma1,sigma2,xi0_in,y);
    % Set in=out to make sure we update:
    xi0_in = xi0_out; 
end

p11
p22
mu1
mu2
sigma1
sigma2
xi0_in

