function [ filteredxi , predictedxi ] = Hamilton_filter(p11,p22,mu,sigma,y)

% Extract length of data
T = length(y);

% Build transition matrix from p11 and p22
P   = [ p11 , 1-p22 ; 1-p11 , p22];


% Run the Hamilton filter
for i=1:T
    % Set the predicted xi
    if i==1
    predictedxi(:,i) = [1;0]; %[ (1-p22) / (2-p11-p22) ; (1-p11)/(2-p11-p22) ];
    else
    predictedxi(:,i) = P * filteredxi(:,i-1);
    end
    likelihood(:,i)  = [ normpdf(y(1,i),mu(1),sigma(1)) ; normpdf(y(1,i),mu(2),sigma(2)) ];
    filteredxi(:,i)  = predictedxi(:,i) .* likelihood(:,i) / ([1,1]*(predictedxi(:,i).*likelihood(:,i)) );
end

% Close the function
end

