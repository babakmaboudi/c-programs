function script()
	load('./data/out.txt')
	out = out';
	dt = 0.0125;	

	X = 50/2000 : 50/2000 : 50;
	for i = 1 : size(out,2)
		plot(X,out(:,i))
		ylim([-1,7])
%		ylim([-1,1])
	
		string = sprintf('%f',dt*i);
		title(string)
		drawnow
	end
	
