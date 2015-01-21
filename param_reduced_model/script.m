function script()
	close all;
	result = load('./data/result.txt');
	result_exact = load('./data/result_exact.txt');

	figure
	hold on
	for i = 1 : length(result)
		plot3(result(i,1),result(i,2),result(i,3),'bx')
		plot3(result_exact(i,1),result_exact(i,2),result_exact(i,3),'r+')
	end
	hold off

	norm(result(end,:) - result_exact(end,:))
