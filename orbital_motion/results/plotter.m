function plotter()

	close all;
	image_expansion_x = load('image_expansion_x.txt');
	image_expansion_y = load('image_expansion_y.txt');
	N=5;
	d=4;

	degrees = pol_degrees(N,d);
	hermit_basis = hermit_pol_coef(N);

	sample = normrnd(0,1,10000,4);
	result_x = [];
	for k = 1 : length(sample)
		u = 0;
		for i = 1 : length(degrees)
			phi = 1;
			for j = 1 : size(degrees,2)
				phi = phi * hermit_pol_eval( hermit_basis , degrees(i,j) , sample(k,j) );
			end
			u = u + image_expansion_x(i) * phi; 
		end
		result_x = [ result_x , u ];
	end

	result_y = [];
	for k = 1 : length(sample)
		u = 0;
		for i = 1 : length(degrees)
			phi = 1;
			for j = 1 : size(degrees,2)
				phi = phi * hermit_pol_eval( hermit_basis , degrees(i,j) , sample(k,j) );
			end
			u = u + image_expansion_y(i) * phi; 
		end
		result_y = [ result_y , u ];
	end

	figure
	hold on
	for i = 1 : length(result_x)
		plot(result_x(i),result_y(i),'+');
	end
	hold off
	axis([2.786e7 2.802e7 -3e6 1.5e6]);



