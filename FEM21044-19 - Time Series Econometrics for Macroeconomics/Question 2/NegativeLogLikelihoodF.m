function [output]=NegativeLogLikelihoodF(F,Q,R,y)

H = eye(4);

% Run the Kalman filter
[~,~,predictedxi,predictedP]=KalmanFilterF(F,Q,R,y);

% Collect a row vector of log likelihood per observation
for t=1:size(y,2)
    Sigma = H'*predictedP(:,:,t)*H + R;
    mu    = H'*predictedxi(:,t);
    LL(t) = log( 1 / sqrt(  det(2*pi*Sigma)  ) * exp(-1/2 *(y(:,t)-mu)' * ((Sigma)\(y(:,t)-mu))));
end

% Sum over all observations (you may use a `burn-in period' of a few observations - but check if it matters)
burninperiod = 1;
output = -sum( LL(burninperiod:end) );               

end