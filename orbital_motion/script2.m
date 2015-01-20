function script2()
	close all;
	A = normrnd(0,1,1000000,1);
	figure;
	hist(A,100)

	figure;
	random = load('random.txt');
	hist(random,100)
