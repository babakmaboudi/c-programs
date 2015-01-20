% Gives the coefficients of the first N hermit polynomilas
% column i corresponts to the coefficients of hermit polynomial of degree i-1, [ a0 a1 ... ai 0 ... 0 ]'

function [coef] = legendre_pol_coef(N)
	
	coef = zeros(N+1,N+1);
	
	coef(1,1) = 1;
	coef(2,2) = 1;

	for i = 3 : N + 1
		n = i - 2;
		coef(:,i) = circshift(coef(:,i-1),1) - n * coef(:,i-2);
	end
