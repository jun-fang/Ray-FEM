xs = 0.5;   ys = 0.2;             % point source location
speed = @(x) ( 1+ 0.2*sin(3*pi*(x(:,1)+0.05)).*sin(0.5*pi*x(:,2)));
cmin = 0.8;

a= 1;
plt = 0;                           % plot solution or not
fquadorder = 3; 
omega = 10*pi;

minwavelength = (2*pi*cmin)/omega;
wpml = 4*minwavelength;        % width of PML
sigmaMax = 30/wpml;                % Maximun absorbtion

NPW = 40;
rh = 1/(NPW*5*round(omega/(2*pi*cmin*5)));
fprintf('omega/(2*pi) = %d,   1/rh = %d,   NPW = %d \n',omega/(2*pi), 1/rh, NPW);
[node,elem] = squaremesh([0,1,0,2],rh);
[ru] = Standard_FEM_PML_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);
figure(1);
FJ_showresult(node,elem,real(ru));
% clear node elem;
% save('test10_reference_solution_20_1.mat','rh','ru');
% test10_caustics;

%% NMLA
Nray = 0;
Rest = 1;
pde = [];
pct = 1/3;
data = 'num';
opt = 0;

minwavelength = (2*pi*cmin)/omega;
r = 4*minwavelength;

m = round(1/rh) + 1;
n = round(2/rh) + 1;


[ux,uy] = num_derivative_2(ru,rh,m,n,2);
% a = 0.5;
% b = 0.5;
% mbd = 0;
% mbd = high_r;
% [mnode,melem] = squaremesh([-a-mbd,a+mbd,-b-mbd,b+mbd],h);
ch = 1/80;
[cnode,celem] = squaremesh([0.2,0.6,0.8,1.2],ch);
cN = size(cnode,1);
cnumray = cell(cN,1);
cray = ex_ray(cnode,xs,ys);

fprintf('NMLA time: \n');
tic;
for i = 1:cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    if d0 < eps 
        cnumray{i} = 0;
    elseif d0 <= r
        cnumray{i} = exp(1i*cray(i,:));
    else
        Rest = 2*d0;
        if d0 <= 2*r
            Rest = 2*d0;
        end
        c0 = speed(cnode(i,:));
        
        [ang] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,node,elem,ru,ux,uy,pde,pct,Nray,data,opt,plt);
        if size(ang,2) > 2
            i
            ang
        end
        cnumray{i} = exp(1i*ang);
    end
end
toc;

figure(2);
ray_field(cnumray,cnode,1);