function script()
	close all
	load('data.txt');
	
	figure;
	for i = 1 : size(data,2)/2
		mat = data(:,i:i+1);
		clf
		hold on
			for j = 1 : size(mat,1)
				plot(mat(j,1),mat(j,2),'*');
			end
		hold off
		drawnow
	end
