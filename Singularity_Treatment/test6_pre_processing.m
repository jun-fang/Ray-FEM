
clear;

load('Babich_CGV.mat');

a = 1/2;
Bh = 1/round( 1/(Bh0/15) );
Bx = -a: Bh : a;  By = -a: Bh : a;
[BX0, BY0] = meshgrid(Bx0, By0);
[BX, BY] = meshgrid(Bx, By);


%% refined phase and amplitude
ttao = interp2(BX0,BY0,tao,BX,BY,'spline'); % refined phase
DD1 = interp2(BX0,BY0,D1,BX,BY,'spline'); % refined amplitude
DD2 = interp2(BX0,BY0,D2,BX,BY,'spline'); % refined amplitude


%% gradient of phase and amplitudes
taox = tao2x ./ (2*tao);   taox(71, 71) = 0;
taoy = tao2y ./ (2*tao);   taoy(71, 71) = 0;
ttaox = interp2(BX0,BY0,taox,BX,BY,'spline'); % refined phase
ttaoy = interp2(BX0,BY0,taoy,BX,BY,'spline'); % refined phase


[D1x,D1y] = num_derivative(D1,Bh0,4);
[D2x,D2y] = num_derivative(D2,Bh0,4);
DD1x = interp2(BX0,BY0,D1x,BX,BY,'spline'); % refined amplitude
DD1y = interp2(BX0,BY0,D1y,BX,BY,'spline'); % refined amplitude
DD2x = interp2(BX0,BY0,D2x,BX,BY,'spline'); % refined amplitude
DD2y = interp2(BX0,BY0,D2y,BX,BY,'spline'); % refined amplitude


%% interpolations
% X = BX(:);  Y = BY(:); 
% Btao = scatteredInterpolant(BX(:), Y(:), ttao(:));
% Btao(0.1, -0.1)

% xy = [0.1, -0.021111];
% uint = interpolation2(Bx, By, ttao, xy)

[cray, cT, cTx, cTy] = eikonal_cgv(1, [0.1, -0.2], [0,0], [BX0(:), BY0(:)]);
norm(cT - tao(:), inf)
norm(cTx - taox(:), inf)
norm(cTy - taoy(:), inf)


[ray, T, Tx, Ty] = eikonal_cgv(1, [0.1, -0.2], [0,0], [BX(:), BY(:)]);
norm(T - ttao(:), inf)
norm(Tx - ttaox(:), inf)
norm(Ty - ttaoy(:), inf)


xy = rand(5,2)/2;
uint = interpolation2(Bx, By, ttao, xy);

[ray, T, Tx, Ty] = eikonal_cgv(1, [0.1, -0.2], [0,0], xy);

norm(uint -T, inf)







%% Construct wave field
% 
% G1 = 1i*(sqrt(pi)/2)*besselh(0,1,omega*ttao); 
% G1 = G1.*DD1;
% 
% G2 = 1i*(sqrt(pi)/2)*exp(1i*pi)*(2*ttao/omega);
% G2 = G2.*besselh(1,1,omega*ttao);
% G2 = G2.*DD2;
% 
% ub = G1 + G2;
% 
% 
% %% Construct gradient of wave field







