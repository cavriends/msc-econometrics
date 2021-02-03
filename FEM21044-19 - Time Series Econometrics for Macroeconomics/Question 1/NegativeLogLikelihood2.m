function [output]=NegativeLogLikelihood2(parameter_vector,y)

% Extract the stuff we need from the input arguments
p11    = parameter_vector(1,1);
p22    = parameter_vector(2,1);
mu     = [ parameter_vector(3,1) ; parameter_vector(4,1) ];
sigma  = [ parameter_vector(5,1) ; parameter_vector(6,1) ];
xi0_in = [ parameter_vector(7,1) ; parameter_vector(8,1) ];

mu1 = mu(1,1);
mu2 = mu(2,1);
sigma1 = sigma(1,1);
sigma2 = sigma(2,1);

% Extract the size of the data
T=size(y,2);

% Run the Hamilton filter
[~,predictedxi] = Hamilton_filter2(p11,p22,mu1,mu2,sigma1,sigma2,xi0_in,y);

% Collect a row vector of log likelihood per observation
for t=1:T
LL(t) =  log( predictedxi(1,t) * my_pdf( y(:,t) , mu1 , sigma1^2 ) + predictedxi(2,t) * my_pdf( y(:,t) , mu2 , sigma2^2 ) );
end

% Sum over all observations and put a minus sign in front
output = - sum( LL(1:end) );
end

