function [ smoothedxi,smoothedP,smoothedPcross,xi0_out,P0_out ] = KalmanSmoother( F , Q , R , y )

% Extract length of the data
T = size(y,2);

% Diffuse initialisation
xi0    = zeros(4,1);
P0     = 10^6*eye(4); % make sure this is the same as the one used in the filter!

% Run the Kalman filter
[xi,P,predictedxi,predictedP]=KalmanFilterF(F,Q,R,y);

% Initialise the smoother at the end
smoothedxi(:,T)   = xi(:,T);
smoothedP(:,:,T)    = P(:,:,T);
        
% Run the Kalman smoother 
for j=1:(T-1)
    % Loop backwards through time
    t = T-j;
    % Run the smoother
    smoothedxi(:,t)   = xi(:,t) + P(:,:,t) * F' * inv( predictedP(:,:,t+1) ) * ( smoothedxi(:,t+1) - predictedxi(:,t+1) ) ;
    smoothedP(:,:,t)    = P(:,:,t)  - P(:,:,t) * F' * inv( predictedP(:,:,t+1) ) * ( predictedP(:,:,t+1) - smoothedP(:,:,t+1) ) * inv( predictedP(:,:,t+1) ) * F * P(:,:,t) ;

    smoothedP(:,:,t) = (smoothedP(:,:,t) + smoothedP(:,:,t)')/2;
end

% Calculate the xi0 and P0 separately
xi0_out = xi0 + P0 * F' * inv( predictedP(:,:,1) ) * ( smoothedxi(:,1) - predictedxi(:,1) ) ;
P0_out  = P0 -  P0 * F' * inv( predictedP(:,:,1) ) * ( predictedP(:,:,1) - smoothedP(:,:,1) ) * inv( predictedP(:,:,1) ) * F * P0 ;

% Run the loop to get the cross terms where smoothedPcross(t)=P_{t,t-1|T}
for t=1:T
    if t==1
        smoothedPcross(:,:,t) = smoothedP(:,:,t) * inv( predictedP(:,:,t) ) * F * P0;
    else
        smoothedPcross(:,:,t) = smoothedP(:,:,t) * inv( predictedP(:,:,t) ) * F * P(:,:,t-1);
    end
end

% Close the function
end

