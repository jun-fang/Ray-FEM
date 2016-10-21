function [] = test1_optimization(omega)


exu = @(p,omega) ( sqrt(omega)*besselh(0,1,omega*sqrt((p(:,1)-2).^2 + (p(:,2)-2).^2)) );

xs = 2;   ys = 2;
pde = Helmholtz_data1;
% global omega;
global a;
   
a= 1;                 
plt = 0;                              % plot solution or not
fquadorder = 3; 
% omega = 200*pi;

NPW = 8;
rh = 1/(NPW*round(omega/(2*pi)));

[rnode,relem] = squaremesh([-a,a,-a,a],rh);
ru = exu(rnode,omega);


%% NMLA 
Rest = 1;
% pde = [];
pct = 1/2; 
data = 'num';

m = size(rnode,1);
m = round(sqrt(m));
n = m;

[ux,uy] = num_derivative_2(ru,rh,m,n,2);
a = 0.5;
b = 0.5;
mbd = 0;
% mbd = high_r;
% [mnode,melem] = squaremesh([-a-mbd,a+mbd,-b-mbd,b+mbd],h);
ch = 1/80;
[cnode,celem] = squaremesh([-a-mbd,a+mbd,-b-mbd,b+mbd],ch);
cN = size(cnode,1);
cnumray1 = zeros(cN,1);
cnumray2 = zeros(cN,1);
cnumray3 = zeros(cN,1);
cnumray4 = zeros(cN,1);


h = rh;
[node,elem] = squaremesh([-a,a,-a,a],h);
exray = pde.ray(node);



Nray = 1;
opt = 1;
fprintf('NMLA time: \n');
tic;
for i = 1:cN
%     i
%     cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    Rest = d0;
    c0 = 1;
    [oangs1,oBs1] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,rnode,relem,ru,ux,uy,pde,pct,Nray,data,0,plt);
    [oangs2,oBs2] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,rnode,relem,ru,ux,uy,pde,pct,Nray,data,opt,plt);
    [angs1,Bs1] = gradient_descent(rnode,ru,omega,x0,y0,c0,oangs1,oBs1);
    [angs2,Bs2] = gradient_descent(rnode,ru,omega,x0,y0,c0,oangs2,oBs2);

    cnumray1(i) = exp(1i*oangs1);
    cnumray2(i) = exp(1i*oangs2);
    cnumray3(i) = exp(1i*angs1);
    cnumray4(i) = exp(1i*angs2);
end
toc;

oangs1
% figure(1);
% ray_field(cnumray1,cnode,1);

numray1 = interpolation(cnode,celem,node,cnumray1);
numray1 = ray_convert(numray1,2);

numray2 = interpolation(cnode,celem,node,cnumray2);
numray2 = ray_convert(numray2,2);

numray3 = interpolation(cnode,celem,node,cnumray3);
numray3 = ray_convert(numray3,2);

numray4 = interpolation(cnode,celem,node,cnumray4);
numray4 = ray_convert(numray4,2);


rayerror1 = norm(exray - numray1);
rayerror2 = norm(exray - numray2);
rayerror3 = norm(exray - numray3)
rayerror4 = norm(exray - numray4);


% figure(2);
% ray_field(numray1,node,2);



%% Ray-FEM
[~,~,~,~,rel_L2_err1] = Ray_FEM_IBC_1(node,elem,omega,pde,numray1,fquadorder,plt);

[~,~,~,~,rel_L2_err2] = Ray_FEM_IBC_1(node,elem,omega,pde,numray2,fquadorder,plt);

[~,~,~,~,rel_L2_err3] = Ray_FEM_IBC_1(node,elem,omega,pde,numray3,fquadorder,plt);

[~,~,~,~,rel_L2_err4] = Ray_FEM_IBC_1(node,elem,omega,pde,numray4,fquadorder,plt);

[~,~,~,~,rel_L2_err_ex] = Ray_FEM_IBC_1(node,elem,omega,pde,exray,fquadorder,plt);


%% print out
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('omega: %d', omega);
fprintf(['\n' '-'*ones(1,80) '\n']);

fprintf('Ray Error: \n');
fprintf('NMLA (basic): %d\n',rayerror1);
fprintf('NMLA (second correction): %d\n',rayerror2);
fprintf('NMLA (basic with optimization): %d\n',rayerror3);
fprintf('NMLA (second correction with optimization): %d\n',rayerror4);


fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Ray-FEM solution relative L2 error: \n');
fprintf('NMLA (basic): %d\n',rel_L2_err1);
fprintf('NMLA (second correction): %d\n',rel_L2_err2);
fprintf('NMLA (basic with optimization): %d\n',rel_L2_err3);
fprintf('NMLA (second correction with optimization): %d\n\n',rel_L2_err4);

fprintf('Exact ray: %d\n',rel_L2_err_ex);
