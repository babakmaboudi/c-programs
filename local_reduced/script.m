function script()
	close all;
	data = load('./data/output.txt');
	exact = load('./data/output2.txt');

	mat = (data(2:end,1:3) - exact(1:end-1,:))';
	sqrt( sum( mat.^2 ) )'

	figure
	hold on

	plot3(data(:,1),data(:,2),data(:,3),'b+')
	plot3(exact(:,1),exact(:,2),exact(:,3),'rx')


	hold off
