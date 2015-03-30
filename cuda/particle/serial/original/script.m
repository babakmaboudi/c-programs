function script()
	close all
	points = load('./data/positions.txt')

	figure;
	hold on
	for i = 1 : length(points)
		plot(points(i,1),points(i,2),'b*');
	end
	hold off
