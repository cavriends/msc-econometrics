function [ xi,P,predictedxi,predictedP ] = KalmanFilterF(F,Q,R,y )
% This function runs the Kalman filter for the scalar AR(1) model plus
% noise with a diffuse prior (roughly)

% Extract lenght of the data
T = size(y,2);
H = eye(4);

for t=1:T
    % Diffuse initialisation
    if t==1
    % Prediction for t=1
    predictedxi(:,t)  = zeros(4,1);
    predictedP(:,:,t)   = 10^6*eye(4);
    % Prediction for any other t
    else
    predictedxi(:,t)  = F*xi(:,t-1);
    predictedP(:,:,t) = F*P(:,:,t-1)*F' + Q;
    
    predictedP(:,:,t) = (predictedP(:,:,t) + predictedP(:,:,t)')/2;
    
    end
    % Update for any t
    xi(:,t)  = predictedxi(:,t) + predictedP(:,:,t)*H*((H'*predictedP(:,:,t)*H + R)\(y(:,t) - H'*predictedxi(:,t)));
    P(:,:,t) = predictedP(:,:,t) - predictedP(:,:,t)*H*((H'*predictedP(:,:,t)*H + R)\H'*predictedP(:,:,t));
    
    P(:,:,t) + (P(:,:,t) + P(:,:,t)')/2;
    % Close the loop over time
end
% Close the function
end
