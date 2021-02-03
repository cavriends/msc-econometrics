function [ xi,P,predictedxi,predictedP ] = KalmanFilterB(parameter_vector,y)
% This function runs the Kalman filter for the scalar AR(1) model plus
% noise with a diffuse prior (roughly)

% Extract lenght of the data
T = size(y,2);

f = [parameter_vector(1,1), parameter_vector(2,1), parameter_vector(3,1), parameter_vector(4,1), parameter_vector(5,1)];
h = [parameter_vector(6,1), parameter_vector(7,1), parameter_vector(8,1), parameter_vector(9,1)];

% Constrain the variances to non-negative values (same as in the
% likelihood)
q = abs([parameter_vector(10,1), parameter_vector(11,1), parameter_vector(12,1), parameter_vector(13,1)]);

% Extract the stuff we need from the input arguments
H = [h;eye(4)];
F = diag(f);
Q = diag([1 q]);

% The Kalman filter for AR(1)+noise
for t=1:T
    % Diffuse initialisation
    if t==1
    % Prediction for t=1
    predictedxi(:,t)  = zeros(5,1);
    predictedP(:,:,t)   = 10^6*eye(5);
    % Prediction for any other t
    else
    predictedxi(:,t)  = F*xi(:,t-1);
    predictedP(:,:,t)   = F*P(:,:,t-1)*F' + Q;
    end
    % Update for any t
    xi(:,t)  = predictedxi(:,t) + predictedP(:,:,t)*H*inv(H'*predictedP(:,:,t)*H)*(y(:,t) - H'*predictedxi(:,t));
    P(:,:,t) = predictedP(:,:,t) - predictedP(:,:,t)*H*inv(H'*predictedP(:,:,t)*H)*H'*predictedP(:,:,t);
    % Close the loop over time
end
% Close the function
end

