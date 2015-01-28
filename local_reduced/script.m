function script()
	close all;
	
	data = load('./data/data.txt');
	centers = load('./data/centers.txt');
	labels = load('./data/labels.txt');
	
	color_code = [ 'b' 'r' 'g' 'k' 'm' 'y' , 'c'];

	figure;
	hold on
	labels = labels + 1;
	for i = 9 : 9
		[ IDX c ] = find(labels == i);
		for j = 1 : length(IDX)
			plot3( data(IDX(j) , 1) , data(IDX(j) , 2) , data(IDX(j) , 3) , [ color_code( mod(i,7) + 1 ), '*'] );
		end
		plot3( centers(i,1) , centers(i,2) , centers(i,3), 'b*' );
	end

	for i = 4 : 4
		[ IDX c ] = find(labels == i);
		for j = 1 : length(IDX)
			plot3( data(IDX(j) , 1) , data(IDX(j) , 2) , data(IDX(j) , 3) , [ color_code( mod(i,7) + 1 ), '*'] );
		end
		plot3( centers(i,1) , centers(i,2) , centers(i,3) , 'b*' );
	end
	for i = 6 : 6
		[ IDX c ] = find(labels == i);
		for j = 1 : length(IDX)
			plot3( data(IDX(j) , 1) , data(IDX(j) , 2) , data(IDX(j) , 3) , [ color_code( mod(i,7) + 1 ), '*'] );
		end
		plot3( centers(i,1) , centers(i,2) , centers(i,3) , 'b*' );
	end
	hold off
