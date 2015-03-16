function script()
	close all;
	data = load('./data/output.txt');
%	exact = load('./data/output2.txt');

	figure
	hold on
	for i = 1 : length(data)
		plot3(data(i,1) , data(i,2) , data(i,3) , 'b+');
%		plot3(exact(i,1) , exact(i,2) , exact(i,3) , 'rx');
	end
	hold off
