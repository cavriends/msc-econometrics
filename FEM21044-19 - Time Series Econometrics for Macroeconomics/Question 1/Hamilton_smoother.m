function [ smoothedxi ] = Hamilton_smoother(p11,p22,mu,sigma,y)

% Extract length of data
T = length(y);

% Build transition matrix from p11 and p22
P   = [ p11 , 1-p22 ; 1-p11 , p22];

% Run the Hamilton filter
[ filteredxi , predictedxi ] = Hamilton_filter(p11,p22,mu,sigma,y);

% Initialise the smoother at time T
smoothedxi(:,T) = filteredxi(:,T);

% Run the `Kim smoother' backwards
for i=1:(T-1)
    t = T - i;
    smoothedxi(:,t) = filteredxi(:,t) .* ( P' * ( smoothedxi(:,t+1) ./ predictedxi(:,t+1) ) );
end
% Notice that in Matlab, .* is elementwise multiplication and ./ is
% elementwise division. Normal matrix multiplication is given by *

% Close the function
end

