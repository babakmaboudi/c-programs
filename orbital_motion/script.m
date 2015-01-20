function script()
	close all;
	monte_carlo = load('monte_carlo.txt');
	
	figure;
	hold on
	for i = 1 : length(monte_carlo)
		plot(monte_carlo(i,1) , monte_carlo(i,2),'+');
	end
	hold off

	mean(monte_carlo(:,1))
	mean(monte_carlo(:,2))	
