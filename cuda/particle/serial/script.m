function script()
	close all
	all = load('./data/positions.txt');

	h = figure;
	for i = 1 : 499
		st = (i-1)*100 + 1;
		en = st + 100;
		points = all(st:en,:);
		hold on
		for j = 1 : 100
			plot(points(j,1),points(j,2),'bo');
		end
		drawnow;
		hold off
		plot(0,0);
	end
	
