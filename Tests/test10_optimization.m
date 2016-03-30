
%% Solution data on a large domain
xs = -0.2;   ys = -0.3;             % point source location
speed = @(x) ( 1+ 0.5*sin(2*pi*x(:,1)));
cmin = 1/2; 
     
a= 1;                 
plt = 0;                              % plot solution or not
fquadorder = 3; 
omega = 60*pi;

wpml = 10*(2*pi*cmin)/omega;        % width of PML
sigmaMax = 50/wpml;                % Maximun absorbtion

NPW = 15;
rh = 1/(NPW*round(omega/(2*pi*cmin)));
fprintf('omega/(2*pi) = %d,   1/h = %d,   NPW = %d \n',omega/(2*pi), 1/rh, NPW);
[rnode,relem] = squaremesh([-a,a,-a,a],rh);
load('test10_reference_solution_20_1.mat');


%% NMLA 
Nray = 0;   
Rest = 1;
pde = [];
pct = 1/2; 
data = 'num';
opt = 0;

minwavelength = (2*pi*cmin)/omega;
r = 4*minwavelength;

m = size(rnode,1);
m = round(sqrt(m));
n = m;

[ux,uy] = num_derivative_2(ru,rh,m,n,2);
a = 0.5;
b = 0.5;
mbd = 0;
% mbd = high_r;
% [mnode,melem] = squaremesh([-a-mbd,a+mbd,-b-mbd,b+mbd],h);
ch = 1/60;
[cnode,celem] = squaremesh([-a-mbd,a+mbd,-b-mbd,b+mbd],ch);

NPW = 7.5;
h = 1/(NPW*round(omega/(2*pi*cmin)));
1/h
[node,elem] = squaremesh([-0.5,0.5,-0.5,0.5],h);
N = size(node,1);
numray = cell(N,1);

cN = size(cnode,1);
cnumray = zeros(cN,1);
cray = ex_ray(cnode,xs,ys,1);
ray = ex_ray(node,xs,ys,1);
cres = zeros(cN,1);
citer = zeros(cN,1);
res1 = zeros(N,1);


Nray = 1;
opt = 1;
fprintf('NMLA time: \n');
tic;
for i = 1:cN
    i
    cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    Rest = d0;
    c0 = speed(cnode(i,:));
    if x0 >= 0 || y0 <= 0.15 % x0 + y0 + 0.1 < 0 || y0 < x0 +  0.5
        [oangs,oBs] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,rnode,relem,ru,ux,uy,pde,pct,Nray,data,opt,plt);
        [angs,Bs,cres(i),citer(i)] = gradient_descent(rnode,ru,omega,x0,y0,c0,oangs,oBs);
    end
    if d0 <= r
        cres(i) = 0;
    end
    cnumray(i) = exp(1i*angs);
end
toc;
numray1 = interpolation(cnode,celem,node,cnumray);
res1 = interpolation(cnode,celem,node,cres);

numray1 = ray_convert(numray1,2);

figure(1);
% ray_field(numray1,node,6);
FJ_showresult(cnode,celem,cres);
 

res = zeros(N,1);   
opt = 0;
pct = 1/3;
tic;
for i = 1:N
    x0 = node(i,1);  y0 = node(i,2);
    d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    numray{i} = numray1(i);
    res(i) = res1(i);
    if d0 <= r
        numray{i} = ray(i);
    elseif x0 < 0 && y0 > 0.15  % (x0 + y0 + 0.1 >= 0) && (y0 >= x0 + 0.5)
        Rest = d0;    
        if d0 <= 2*r
            Rest = 2*d0;
        end   
        c0 = speed(node(i,:));
        Nray = 0;
        [oangs,oBs] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,rnode,relem,ru,ux,uy,pde,pct,Nray,data,opt,plt);
        if and(and(x0 + y0 + 0.1 >= 0, y0 >= x0 + 0.5), (size(oangs) == 1))
            pct = 1/5;
            [oangs,oBs] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,rnode,relem,ru,ux,uy,pde,pct,Nray,data,opt,plt);
            pct = 1/3;
        end
%         if sum(and(oangs/pi>0.75, oangs/pi<1)) == 0
%             oangs = [oangs, 0.8536*pi];
%             oBs = [oBs, 0.0060 + 0.0104i]; 
%         end
        [angs,Bs,res(i),~] = gradient_descent(rnode,ru,omega,x0,y0,c0,oangs,oBs);
        if size(angs,2) > 20
            temp1 = abs(Bs);
            temp2 = sort(temp1);
            ind = find(temp2>temp2(end-1)-10*eps);
            angs = angs(ind);
            Bs = Bs(ind);
            [angs2,Bs2,res2,~] = gradient_descent(rnode,ru,omega,x0,y0,c0,oangs,oBs);
            if res2 < res(i)
                angs = angs2;
            end
        end
        
        numray{i} = exp(1i*angs);
    end
end
toc;

figure(2);
ray_field(numray,node,6);
figure(3);
FJ_showresult(node,elem,res);    



%% Ray_FEM
wpml = 10*(2*pi*cmin)/omega;        % width of PML
sigmaMax = 50/wpml; 
[u2] = Ray_FEM_PML_22_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,numray,fquadorder,plt);
% figure(4);
% FJ_showresult(node,elem,real(u2));


% [us] = Standard_FEM_PML_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);



%% Reference
fh = 1/2000;
[fnode,felem] = squaremesh([-0.5,0.5,-0.5,0.5],fh);
[uf] = Standard_FEM_PML_PointSource(fnode,felem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);
% figure(5);
% FJ_showresult(fnode,felem,real(uf));

%% plot
a = 0.5;   
[X,Y] = meshgrid(-a:h:a,-b:h:b);
[m,n] = size(X);

u2 = reshape(u2,m,n);
% us = reshape(us,m,n);


fn = round(1/fh) + 1;
uf = reshape(uf,fn,fn);


% figure(7);
% hold off;
% dy = 0.8;
% yn = round(dy/h) + 1;
% xx = X(yn,:);
% uu = u2(yn,:);
% ryn = round(dy/fh) + 1;
% rxx = -0.5:fh:0.5;
% ruu = uf(ryn,:);
% plot(xx,real(uu),'ro-');
% hold on;
% ss = us(yn,:);
% plot(xx,real(ss),'bs-');
% hold on;
% plot(rxx,real(ruu),'k');
% xlabel('x');
% ylabel('Wavefield');
% legend('Ray-FEM solution','standard FEM solution','Reference solution','LOCATION','Best');

figure(9);
hold off;
dy = 0.8;
yn = round(dy/h) + 1;
xx = X(yn,:);
uu = u2(yn,:);
ryn = round(dy/fh) + 1;
rxx = -0.5:fh:0.5;
ruu = uf(ryn,:);
plot(xx,real(uu),'ro-');
hold on;
plot(rxx,real(ruu),'k');
xlabel('x');
ylabel('Wavefield');
legend('Ray-FEM solution','Reference solution','LOCATION','Best');



figure(10);   
hold off;
dy = 0.85;
yn = round(dy/h) + 1;
xx = X(yn,:);
uu = u2(yn,:);
ryn = round(dy/fh) + 1;
rxx = -0.5:fh:0.5;
ruu = uf(ryn,:);
plot(xx,real(uu),'ro-');
hold on;
plot(rxx,real(ruu),'k');
xlabel('x');
ylabel('Wavefield');
legend('Ray-FEM solution','Reference solution','LOCATION','Best');




