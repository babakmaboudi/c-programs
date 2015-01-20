function script()
	close all
	mat = load('./data/result.txt');
	exact = load('./data/result_exact.txt')

	figure
	hold on
	for i = 1 : length(mat)
		plot3(mat(i,1),mat(i,2),mat(i,3),'b+')
		plot3(exact(i,1),exact(i,2),exact(i,3),'rx')
	end
	hold off

	norm( mat(end,:) - exact(end,:) )
       whos	
