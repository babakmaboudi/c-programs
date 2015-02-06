function script()
	close all;
	data = load('./data/output.txt');

	figure
	hold on
	for i = 1 : length(data)
		plot3(data(i,1) , data(i,2) , data(i,3) , 'b*');
	end
	hold off
