function [output]=NegativeLogLikelihoodA(parameter_vector,y)

% Extract the stuff we need from the input arguments
F = parameter_vector(1,1);
Q = parameter_vector(2,1);
R = parameter_vector(3,1);

% Run the Kalman filter
[~,~,predictedxi,predictedP]=KalmanFilterA(parameter_vector,y);

% Collect a row vector of log likelihood per observation
for t=1:size(y,2)
    Sigma(t) = predictedP(t)+R;
    mu(t)    = predictedxi(t);
    LL(t)    = log( 1 / sqrt(  det(2*pi*Sigma(t))  ) * exp( -1/2 * (y(t)-mu(t))' * inv(Sigma(t)) *  (y(t)-mu(t))  )   );
end

% Sum over all observations (you may use a `burn-in period' of a few observations - but check if it matters)
burninperiod = 1;
output = - sum( LL(burninperiod:end) )        ;               

end

