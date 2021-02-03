function [ output ] = my_pdf( y , mu , Sigma )

output = 1/sqrt( det(2 * pi * Sigma)) * exp( -1/2 * (y-mu)' * inv(Sigma) * (y-mu) );

end

