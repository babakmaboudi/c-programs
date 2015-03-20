function script2()
	A = load('./data/a.txt');
	[nrows , ncols] = size(A);

	j2k = [ zeros( ncols/2 , ncols/2 ) , eye( ncols/2 ) ; - eye( ncols/2 ) , zeros( ncols/2 , ncols/2 ) ];
	j2n = [ zeros( nrows/2 , nrows/2 ) , eye( nrows/2 ) ; - eye( nrows/2 ) , zeros( nrows/2 , nrows/2 ) ];

	A_inv = j2k'*A'*j2n;
	L = zeros(6,6);
	L(1,4) = 1;
	L(2,5) = 1;
	L(3,6) = 1;
	A_inv * L * A

