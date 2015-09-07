function str = give_name(index)
	counter = 0;

	dummy = index;
	while( fix(dummy/10) > 0 )
		dummy = fix(dummy/10)
		counter = counter +1;
	end

	str = 'imag';
	for i = 1 : 5 - counter
		str = [ str , '0' ];
	end
	temp = int2str(index);

	str = [str , temp];
	str = [str , '.png'];
