function script()
	A = load('./data/solution_reduced.txt');
	for i=1 : size(A,2)
		plot(A(:,i));
		ylim([-1,1])
		pause(0.01);
	end

