function draws = inv_gaussian_generator( mu, lambda, p  )
% this function is copied from Luc Devroye  pag 149.
%================================================================
%mu is the mu parameter of the distribution
%lambda is the lambda parameter of the distribution 
% p is the dimension of mu and lambda and determines the size of the draws.
% we read each component of mu and lambda and return a draw for each mu and
% lambda.
% to draw many values from the same distribution (good to check code) input
% a large vector with copies of the same mu and lambda.

draws = zeros(p,1);

for i = 1: p
    
    v = randn(1,1);
    y = v^2;    
    m = mu( i );
    l = lambda( i );
        
    x = m + ( ( (m^2.0) * y ) / ( 2.0 * l) ) -...
          ( ( m /(2.0 * l)) * ( sqrt( ( 4.0 * m * l * y) + ( (m^2.0) * (y^2.0 )) ) ) ) ;
    
    u = rand(1, 1);
    
    if   ( u <= ( (m/(m + x) ) ) ) 
        
        keep = x;
       
    else
        
         keep = ( (m^2.0) /x ); 
        
    end
    
    draws( i , 1 ) = keep;
    
end
   
end

% theoretical first moment is mu and the second is mu^3/lambda.
% thus the mean and variance of random draws should be close to those
% theoretical values ;o)