function script()
	close all;

	data = load('./data/out.txt');

	for i = 1 : length(data)
		instance = data(i,:);
		plot(instance)
		ylim([0,1])
		pause(0.01)
	end
