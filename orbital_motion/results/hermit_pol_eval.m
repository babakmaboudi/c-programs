% Computes the value of legendre polynomial of degree n at 
% point x with the givven set of coefficients coef
function value = hermit_pol_eval(coef,n,x)
	value = zeros( 1 , length(x) );
	p = ones( 1 , length(x) );
	for i = 1 : n + 1
		value = value + coef( i , n + 1 ) .* p;
		p = p .* x;
	end
