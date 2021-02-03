function [output]=NegativeLogLikelihood(parameter_vector,y)

% Extract the stuff we need from the input arguments
p11   = parameter_vector(1,1);
p22   = parameter_vector(2,1);
mu    = [ parameter_vector(3,1) ; parameter_vector(4,1) ];

% Added the constraint of non-negative variances
sigma = abs([ parameter_vector(5,1) ; parameter_vector(6,1) ]);

% Run the Hamilton filter
[~,predictedxi] = Hamilton_filter(p11,p22,mu,sigma,y);

% Collect a row vector of log likelihood per observation
LL =  log( predictedxi(1,:) .* normpdf( y , mu(1) , sigma(1) ) + predictedxi(2,:) .* normpdf( y , mu(2) , sigma(2) ) );

% Sum over all observations (or you may use a burn-in period)
burninperiod=1;
output = - sum( LL(burninperiod:end) );               

end

