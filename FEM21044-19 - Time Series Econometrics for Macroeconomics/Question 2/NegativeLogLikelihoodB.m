function [output]=NegativeLogLikelihoodB(parameter_vector,y)

% Extract the stuff we need from the input arguments
f = [parameter_vector(1,1), parameter_vector(2,1), parameter_vector(3,1), parameter_vector(4,1), parameter_vector(5,1)];
h = [parameter_vector(6,1), parameter_vector(7,1), parameter_vector(8,1), parameter_vector(9,1)];

% Constrain the variances to be non-negative values (same as in the filter)
q = abs([parameter_vector(10,1), parameter_vector(11,1), parameter_vector(12,1), parameter_vector(13,1)]);

% Extract the stuff we need from the input arguments
H = [h;eye(4)];
F = diag(f);
Q = diag([1 q]);

% Run the Kalman filter
[~,~,predictedxi,predictedP]=KalmanFilterB(parameter_vector,y);

% Collect a row vector of log likelihood per observation
for t=1:size(y,2)
    Sigma = H'*predictedP(:,:,t)*H;
    mu    = H'*predictedxi(:,t);
    LL(t) = log( 1 / sqrt(  det(2*pi*Sigma)  ) * exp(-1/2 *(y(:,t)-mu)' * ((Sigma)\(y(:,t)-mu))));
end

% Sum over all observations (you may use a `burn-in period' of a few observations - but check if it matters)
burninperiod = 1;
output = -sum( LL(burninperiod:end) );               

end