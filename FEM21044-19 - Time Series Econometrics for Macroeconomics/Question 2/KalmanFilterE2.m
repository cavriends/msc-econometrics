function [ xi,P,predictedxi,predictedP ] = KalmanFilterE2(parameter_vector,y)
% This function runs the Kalman filter for the scalar AR(1) model plus
% noise with a diffuse prior 

% Extract lenght of the data
T = size(y,2);

f = [parameter_vector(1,1), parameter_vector(2,1), parameter_vector(3,1), parameter_vector(4,1), parameter_vector(5,1)];
h = [parameter_vector(6,1), parameter_vector(7,1), parameter_vector(8,1), parameter_vector(9,1)];

% Constrain the variances to non-negative values (same as in the
% likelihood)
q = abs([parameter_vector(10,1), parameter_vector(11,1), parameter_vector(12,1), parameter_vector(13,1)]);
r = abs([parameter_vector(14,1), parameter_vector(15,1), parameter_vector(16,1), parameter_vector(17,1)]);

% Extract the stuff we need from the input arguments
H = [h;eye(4)];
F = diag(f);
Q = diag([1 q]);
R = diag(r);

% Diffuse initialisation
xi0  = zeros(5,1);
P0   = 10^6*eye(5); % make sure this is the same as the one used in the smoother!

% The Kalman filter
for t=1:T
    if t<2
    % The first prediction is based on xi0 and P0
    predictedxi(:,t)    = F * xi0;
    predictedP(:,:,t)     = F * P0 * F'+ Q;
    else
    % Further predictions are based on previous updates
    predictedxi(:,t)    = F * xi(:,t-1);
    predictedP(:,:,t)     = F * P(:,:,t-1) * F' + Q;
    end
    % Update
    xi(:,t)  = predictedxi(:,t)  + predictedP(:,:,t)* H * inv(H'*predictedP(:,:,t)*H + R)*( y(:,t) - H'*predictedxi(:,t));
    P(:,:,t)   = predictedP(:,:,t)   - predictedP(:,:,t) * H * inv( H'*predictedP(:,:,t)*H + R)* H' * predictedP(:,:,t);
    % Close the loop over time
    end
% Close the function  
end

