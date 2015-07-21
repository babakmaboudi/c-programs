function script()
	close all
	all = load('./data/positions.txt');
	vel = load('./data/velocity.txt');

	h = figure;
	for i = 1 : 999
		st = (i-1)*100 + 1;
		en = st + 99;
		points = all(st:en,:);
		stn = i*100 + 1;
		enn = stn + 99;
		pointsn = all(stn:enn,:);
		hold on
		rectangle('Position',[-0.003,-0.003,0.223607+0.006,0.223607+0.006])
		for j = 1 : 100
%			viscircles([points(j,1),points(j,2)],0.003)
%			plot(points(j,1),points(j,2),'bo');
%			filledCircle([points(j,1),points(j,2)],0.0035,100,'b');
			vx = [points(j,1), pointsn(j,1)]; 
			vy = [points(j,2), pointsn(j,2)]; 
			arrowline(vx,vy);
		end
		axis equal
		xlim([-0.05 0.223607+0.05])
		ylim([-0.05 0.223607+0.05])
		
		drawnow;
		hold off
		str = give_name(i-1);
		saveas(h,str);
		plot(0,0);
	end


function h = filledCircle(center,r,N,color)
	THETA=linspace(0,2*pi,N);
	RHO=ones(1,N)*r;
	[X,Y] = pol2cart(THETA,RHO);
	X=X+center(1);
	Y=Y+center(2);
	h=fill(X,Y,color);

function str = give_name(index)
	counter = 0;

	dummy = index;
	while( fix(dummy/10) > 0 )
		dummy = fix(dummy/10);
		counter = counter +1;
	end

	str = './out/imag';
	for i = 1 : 5 - counter
		str = [ str , '0' ];
	end
	temp = int2str(index);

	str = [str , temp];
	str = [str , '.png'];
