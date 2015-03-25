function script2()
	mat1r = load('out_real1.txt');
	mat1i = load('out_imag1.txt');
	mat2r = load('out_real2.txt');
	mat2i = load('out_imag2.txt');

	mat1 = mat1r + 1i*mat1i;
	mat2 = mat2r + 1i*mat2i;

	mat1 * mat2
