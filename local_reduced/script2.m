function script2()
	close all
	result = load('./data/result_reduced.txt');
	exact = load('./data/result_exact.txt');

	figure;
	hold on
	for i = 1 : length(result)
		plot3( result(i,1) , result(i,2) , result(i,3) , 'bx' );
		plot3( exact(i,1) , exact(i,2) , exact(i,3) , 'r+' );
	end
	hold off
%
%	result(end,:)
%	exact(end,:)
	norm( result(end,:) - exact(end,:) )
