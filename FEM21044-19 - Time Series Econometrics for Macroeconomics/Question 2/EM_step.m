function [ F_out,Q_out,R_out ] = EM_step( F_in , Q_in , R_in , y )

% When running EM, we are really looking for a fixed point, such that "output = input" holds.
% If this holds, then we are done.... but generally it doesn't, which is why
% we iterate. In theory we should iterate an infinite number of times to
% get convergence. 
% However, in practice, we use some convergence criterion, which says that we stop
% when the ouput and input differ by 'only a little bit'. 
% For example, we could stop when the input and output are the same up to M decimal places, 
% for some pre-specified precision level M.

%% Extract length of the data
T = size(y,2);

% Diffuse initialisation
xi0       = zeros(4,1);
P0        = 10^6*eye(4); % make sure this is the same as in the filter and smoother

% In this demo code, it is assumed that the initial state is
% truly random: anything could happen. Thus each EM step starts with a
% diffuse prior. Since we are comparing with maximum likelihood, which is also initiated
% with a diffuse prior, it is only fair that we do the same for EM. 

%% E-step: Run the Kalman smoother

[smoothedxi,smoothedP,smoothedPcross,xi0_out,P0_out] = KalmanSmoother(F_in,Q_in,R_in,y) ;

% The Kalman smoother has more outputs than usual. It outputs the cross
% covariance smoothedPcross, as well as xi0_out and P0_out. 

%% Optimise over F
for t=1:T
    if t==1
        numerator(:,:,t) = smoothedxi(:,t)* xi0_out'          + smoothedPcross(:,:,t) ;
    else 
        numerator(:,:,t) = smoothedxi(:,t)*smoothedxi(:,t-1)' + smoothedPcross(:,:,t) ; 
    end
end

for t=1:T
    if t==1
        denominator(:,:,t) = xi0_out * xi0_out'                     + P0_out ; 
    else
        denominator(:,:,t) = smoothedxi(:,t-1) * smoothedxi(:,t-1)' + smoothedP(:,:,t-1) ;
    end
end

F_out = sum(numerator,3) / sum(denominator,3);
% note that A / B is the same as A * inv(B)

clear numerator denominator

%% Optimise over Q

for t=1:T
    if t==1
        Q_estimate(:,:,t) = smoothedxi(:,t) * smoothedxi(:,t)' + smoothedP(:,:,t) - F_out * ( xi0_out           * smoothedxi(:,t)'  + smoothedPcross(:,:,t)' ) - ( smoothedxi(:,t) * xi0_out'        + smoothedPcross(:,:,t) ) * F_out' + F_out * ( xi0_out * xi0_out'                 + P0_out         ) * F_out' ;
    else
        Q_estimate(:,:,t) = smoothedxi(:,t) * smoothedxi(:,t)' + smoothedP(:,:,t) - F_out * ( smoothedxi(:,t-1)  * smoothedxi(:,t)'  + smoothedPcross(:,:,t)' ) - ( smoothedxi(:,t) * smoothedxi(:,t-1)' + smoothedPcross(:,:,t) ) * F_out' + F_out * ( smoothedxi(:,t-1) * smoothedxi(:,t-1)' + smoothedP(:,:,t-1) ) * F_out' ;
    end
end

Q_mean = mean(Q_estimate,3) ;

Q_out = (Q_mean + Q_mean')/2;

clear Q_estimate

%% Optimise over R

for t=1:T
R_estimate(:,:,t) =  y(:,t)*y(:,t)' - smoothedxi(:,t) * y(:,t)'  - y(:,t) * smoothedxi(:,t)' + ( smoothedxi(:,t) * smoothedxi(:,t)' + smoothedP(:,:,t) ) ;
end

R_mean = mean(R_estimate,3) ;

R_out = (R_mean + R_mean')/2;

clear R_estimate

%%Close the function
end

