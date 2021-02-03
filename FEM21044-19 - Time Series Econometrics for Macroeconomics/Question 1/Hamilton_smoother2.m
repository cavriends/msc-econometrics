function [ smoothedxi ,xi0_out , Pstar ] = Hamilton_smoother2(p11,p22,mu1,mu2,sigma1,sigma2,xi0_in,y)

% Extract length of data
T = length(y);

% Build transition matrix from p11 and p22
P   = [ p11 , 1-p22 ; 1-p11 , p22];

% Run the Hamilton filter
[ filteredxi , predictedxi ] = Hamilton_filter2(p11,p22,mu1,mu2,sigma1,sigma2,xi0_in,y);

% Initialise the smoother at time T
smoothedxi(:,T) = filteredxi(:,T);

% Run the `Kim smoother' backwards
for i=1:(T-1)
    t = T - i;
    smoothedxi(:,t) = filteredxi(:,t) .* ( P' * ( smoothedxi(:,t+1) ./ predictedxi(:,t+1) ) );
end

% Separately compute the smoothed state at time zero
xi0_out = xi0_in .* ( P' * ( smoothedxi(:,1) ./ predictedxi(:,1) ) );

% Get the cross terms
for t=1:T
    if t==1
    %Pstar(:,:,t) = P.* (xi0_out         * xi0_in') ./ ( xi0_in * [1,1]) ;  
    % correction suggested by your fellow student
    Pstar(:,:,t) = P.* (smoothedxi(:,1) * xi0_in') ./ ( predictedxi(:,1) * [1,1]);
    else
    Pstar(:,:,t) = P.* (smoothedxi(:,t) * filteredxi(:,t-1)') ./ ( predictedxi(:,t) * [1,1]) ;
    end

% Notes:

% Notice that in Matlab, .* is elementwise multiplication and ./ is
% elementwise division. Normal matrix multiplication is given by *

% To avoid numerical inaccuracies, we can enforce the normalisation if we
% want:
smoothedxi = smoothedxi ./ repmat( sum(smoothedxi) , 2, 1);

% To avoid numerical inaccuracies, we may normalise xi0_out in a hard
% fashion:
xi0_out = xi0_out ./ repmat( sum(xi0_out) ,2 ,1);

% To avoid numerical inaccuracies, we can enfore the normalisation
Pstar = Pstar ./ sum(sum(Pstar));

% Close the function
end

