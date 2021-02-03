function [ xi,P,predictedxi,predictedP ] = KalmanFilterA( parameter_vector,y )
% This function runs the Kalman filter for the scalar AR(1) model plus
% noise with a diffuse prior (roughly)

% Extract lenght of the data
T = size(y,2);

% Extract the stuff we need from the input arguments
F = parameter_vector(1,1);
Q = parameter_vector(2,1);
R = parameter_vector(3,1);

% Diffuse initialisation
xi0  = 0;
P0   = 10^6; % make sure this is the same as the one used in the smoother!

% The Kalman filter for AR(1)+noise
for t=1:T
    % Diffuse initialisation
    if t<2
    % Prediction for t=1
    predictedxi(t)  = F * xi0; %0;
    predictedP(t)   = F * P0 * F' + Q; %10^6;
    % Prediction for any other t
    else
    predictedxi(t)  = F * xi(t-1);
    predictedP(t)   = F * P(t-1) * F' + Q;
    end
    % Update for any t
    xi(t)  = predictedxi(t)  + predictedP(t) * 1/( predictedP(t) + R ) * ( y(t) - predictedxi(t) );
    P(t)   = predictedP(t)   - predictedP(t) * 1/( predictedP(t) + R ) * predictedP(t);
    % Close the loop over time
end
% Close the function
end

