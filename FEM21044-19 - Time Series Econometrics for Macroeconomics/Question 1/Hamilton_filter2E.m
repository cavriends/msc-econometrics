function [ filteredxi , predictedxi] = Hamilton_filter2(p11,p22,mu1,mu2,sigma1,sigma2,xi0_in,y)

% Extract length of data
T = size(y,2);

% Build transition matrix from p11 and p22
P   = [ p11 , 1-p22 ; 1-p11 , p22];

% Predict using the user-defined input of the state at time zero, which we
% call xi0_in (zero for time zero, "in" for input by us)
predictedxi(:,1) = P * xi0_in;

for i=1:T
   likelihood(:,i)   = [ my_pdf(y(:,i),mu1,sigma1) ; my_pdf(y(:,i),mu2,sigma2) ];
   filteredxi(:,i)   = predictedxi(:,i) .* likelihood(:,i) / (predictedxi(:,i)' * likelihood(:,i) );
   predictedxi(:,i+1)= P * filteredxi(:,i) ;
end

% Delete the last prediction, because we want filteredxi and
% predictedxi to have the same length
predictedxi = predictedxi(:,1:T);

% To avoid possible numerical inaccuracies, it is possible to make sure
% both xi are normalised
filteredxi  = filteredxi ./ repmat( sum(filteredxi) , 2, 1);
predictedxi = predictedxi ./ repmat( sum(predictedxi) , 2, 1);

% Close the function
end

