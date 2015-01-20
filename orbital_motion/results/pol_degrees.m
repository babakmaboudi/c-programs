function mat = pol_degrees(N,d)
	mat = zeros(1,d);
	for n = 1 : N
		instance = zeros(1,d);
		instance(1) = n;
	
		mat = [mat ; instance];
		ptr = 1;
		while(instance(end) ~= n)
			while(ptr ~= d)
				instance(ptr) = instance(ptr) - 1;
				ptr = ptr + 1;
				instance(ptr) = instance(ptr) + 1;
				mat = [ mat ; instance ];
			end
			temp = instance(ptr);
			instance(ptr) = 0;
			while( (instance(ptr) == 0) && (ptr > 1) )
				ptr = ptr - 1;
			end
			if instance(ptr) == 0
				instance(end) = n;
			else
				instance(ptr) = instance(ptr) - 1;
				ptr = ptr + 1;
				instance(ptr) = temp + 1;
				mat = [ mat ; instance ];
			end
		end
	end
