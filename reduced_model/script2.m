function script2()
	Umat = load('./data/U.txt');
	phi = load('./data/phi.txt');
	init_cond = load('./data/init_cond.txt');
	init_cond_red = load('./data/init_cond_red.txt');

	snap = load('./data/all_result.txt');

	[U,S,V] = svd(snap');
	Phimat = U(:,1:4);

	a = Phimat' * init_cond;
	Phimat * a - phi * ( phi' * init_cond )
