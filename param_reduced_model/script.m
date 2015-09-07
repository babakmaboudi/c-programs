function script()
	close all;
	data = load('./data/result.txt');
	exact = load('./data/result_exact.txt');

	figure
	hold on

	plot3(data(:,1),data(:,2),data(:,3),'b')
	plot3(exact(:,1),exact(:,2),exact(:,3),'r')


	hold off
