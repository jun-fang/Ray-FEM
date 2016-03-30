% function [] = test1_optimization(omega)

pde = Helmholtz_data3;
global omega;
global a;
   
a= 1;                 
plt = 0;                              % plot solution or not
fquadorder = 4; 
omega = 100*pi;

NPW = 8;
rh = 1/(NPW*round(omega/(2*pi)));

[rnode,relem] = squaremesh([-a,a,-a,a],rh);
ru = pde.ex_u(rnode);


%% NMLA 
Rest = 3.4;
pct = 1/8; 
data = 'num';

m = size(rnode,1);
m = round(sqrt(m));
n = m;

[ux,uy] = num_derivative_2(ru,rh,m,n,2);
a = 0.5;
b = 0.5;
mbd = 0;
ch = 1/80;
[cnode,celem] = squaremesh([-a-mbd,a+mbd,-b-mbd,b+mbd],ch);
cN = size(cnode,1);
cnumray1 = zeros(cN,4);   % basic NMLA
cnumray2 = zeros(cN,4);   % basic NMLA with optimization


h = rh;
[node,elem] = squaremesh([-a,a,-a,a],h);
exray = pde.ray(node);



Nray = 0;
opt = 1;
fprintf('NMLA time: \n');
tic;
for i = 1:cN
%     i
%     cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    c0 = 1;
    [oangs1,oBs1] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,rnode,relem,ru,ux,uy,pde,pct,Nray,data,0,plt);
    [angs1,Bs1] = gradient_descent(rnode,ru,omega,x0,y0,c0,oangs1,oBs1);
    if length(oangs1)~=4
        fprintf('Not four ray directions: %d',oangs1);
        break;
    end

    cnumray1(i,:) = exp(1i*oangs1);
    cnumray2(i,:) = exp(1i*angs1);
end
toc;


% figure(1);
% ray_field(cnumray1,cnode,1);

numray1 = interpolation(cnode,celem,node,cnumray1);
numray1 = ray_convert(numray1,2);

numray2 = interpolation(cnode,celem,node,cnumray2);
numray2 = ray_convert(numray2,2);

dray1 = exray - numray1;
dray2 = exray - numray2;


rayerror1 = norm(dray1(:))
rayerror2 = norm(dray2(:))


% figure(2);
% ray_field(numray1,node,2);



%% Ray-FEM
[~,~,~,~,rel_L2_err1] = Ray_FEM_IBC_1(node,elem,omega,pde,numray1,fquadorder,plt);

[~,~,~,~,rel_L2_err2] = Ray_FEM_IBC_1(node,elem,omega,pde,numray2,fquadorder,plt);

[~,~,~,~,rel_L2_err_ex] = Ray_FEM_IBC_1(node,elem,omega,pde,exray,fquadorder,plt);


%% print out
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('omega: %d', omega);
fprintf(['\n' '-'*ones(1,80) '\n']);

fprintf('Ray Error: \n');
fprintf('NMLA (basic): %d\n',rayerror1);
fprintf('NMLA (basic with optimization): %d\n',rayerror2);


fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Ray-FEM solution relative L2 error: \n');
fprintf('NMLA (basic): %d\n',rel_L2_err1);
fprintf('NMLA (basic with optimization): %d\n',rel_L2_err2);

fprintf('Exact ray: %d\n',rel_L2_err_ex);
