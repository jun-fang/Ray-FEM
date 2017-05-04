a = 1:4;
a = 10.^a;

b1 = abs(besselh(0,1,a));
b2 = abs(besselh(1,1,a));
b3 = abs(besselh(2,1,a));

figure(123);
subplot(1,3,1);
showrate(a,b1);
subplot(1,3,2);
showrate(a,b2);
subplot(1,3,3);
showrate(a,b3);