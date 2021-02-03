function [ xi,P,predictedxi,predictedP ] = KalmanFilterF2( F,Q,R,y )
% This function runs the Kalman filter for the scalar AR(1) model plus
% noise with a diffuse prior 

% Extract lenght of the data
T = size(y,2);

% Diffuse initialisation
xi0  = zeros(4,1);
P0   = 10^6*eye(4); % make sure this is the same as the one used in the smoother!

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
    xi(:,t)  = predictedxi(:,t)  + predictedP(:,:,t) * inv( predictedP(:,:,t) + R ) * ( y(:,t) - predictedxi(:,t) );
    P(:,:,t)   = predictedP(:,:,t)   - predictedP(:,:,t) * inv( predictedP(:,:,t) + R ) * predictedP(:,:,t);
    % Close the loop over time
    end
% Close the function  
end

