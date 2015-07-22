function script()
	close all
	load('data.txt');
	
%	figure;
%	hold on
%	for i = 1 : size(data,1)
%		plot3(data(i,1),data(i,2),data(i,3),'*')
%	end
%	hold off

	[a,b] = max(data(:,1:6))

	figure;
	for i = 1 : size(data,2)/3
		mat = data(:,i:i+2);
		clf
		hold on
			for j = 1 : size(mat,1)
				plot3(mat(j,1),mat(j,2),mat(j,3),'*');
			end
		hold off
		view([1,1,1])
		drawnow
	end
