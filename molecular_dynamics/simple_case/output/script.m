function script()
	load('data.txt')
	data = data';

	[U S V] = svd(data);

	diag(S)
