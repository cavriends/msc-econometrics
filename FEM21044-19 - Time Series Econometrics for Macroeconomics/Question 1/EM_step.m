function [p11_out,p22_out,mu1_out,mu2_out,sigma1_out,sigma2_out,xi0_out]=EM_step(p11,p22,mu1,mu2,sigma1,sigma2,xi0_in,y)

%
T=size(y,2);

% Expectation step: run the smoother
[smoothedxi,xi0_out,Pstar] = Hamilton_smoother2(p11,p22,mu1,mu2,sigma1,sigma2,xi0_in,y);
    
% Extract the star quantities
for t=1:T
    p11star(t) = Pstar(1,1,t);
    p12star(t) = Pstar(1,2,t);
    p1star(t) =  p11star(t) + p12star(t);
    p21star(t) = Pstar(2,1,t);
    p22star(t) = Pstar(2,2,t);
    p2star(t)  = p21star(t) + p22star(t);
end

% Maximisation step, analytic solution for the parameters:
p11_out   = sum( p11star ) / (xi0_out(1) + sum( p1star(1:end-1) )) ;
p22_out   = sum( p22star ) / (xi0_out(2) + sum( p2star(1:end-1) )) ;
mu1_out   = sum( p1star .* y ) / sum( p1star ); 
mu2_out   = sum( p2star .* y ) / sum( p2star ); 
sigma1_out= sqrt( sum( p1star .* ( y - mu1 ).^2 ) / sum(p1star) );
sigma2_out= sqrt( sum( p2star .* ( y - mu2 ).^2 ) / sum(p2star) );

% Close the function
end

