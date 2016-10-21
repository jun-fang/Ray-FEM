%% One point source inside domain (heterogeneous medium)

xc = 0;   yc = 0;           % point source location
% speed = @(x) (3 - 2.5*exp( -((x(:,1)+1/8).^2 + (x(:,2)-0.1).^2)/0.8^2 )) ;    % medium speed
    
   
plt = 0;                         % plot solution or not
fquadorder = 3;                  % numerical quadrature order
Nray = 1;                        % no ray crosssing, number of ray direction is 1
        
high_omega = 100*pi;              % high frequency
low_omega = 2*sqrt(high_omega);    % low frequency
NPW = 16;                          % number of grid points per wavelength

h = 1/2^(round(log2((NPW*high_omega)/(2*pi))));    %% 1/h should be a mutiple of 40!!
ch = 4*h;% 1/2^(round(log2((NPW*low_omega)/(2*pi))));

high_wl = 2*pi/high_omega;
high_wpml = ch*round(high_wl/ch);     % width of PML
high_sigmaMax = 25/high_wpml;                  % Maximun absorbtion
low_wl = 2*pi/low_omega;
low_wpml = low_wl;     % width of PML
low_sigmaMax = 25/low_wpml;                  % Maximun absorbtion
dr = 1.2*low_wl;

% speed = @(x) (1 - 0.5*exp( -((x(:,1)-xc).^2 + (x(:,2)-yc).^2)/0.8^2 ))...
%     .*(sqrt((x(:,1)-xc).^2 + (x(:,2)-yc).^2)>dr) +...
%     (1 - 0.5*exp( -(dr^2)/0.8^2 )).*(sqrt((x(:,1)-xc).^2 + (x(:,2)-yc).^2)<=dr);    % medium speed

speed = @(x) 1./sqrt((x(:,1)+20).^2 + (x(:,2)+10).^2); 

sm_a = 1/2;
md_a = 1/2 + ch*ceil(high_wpml/ch);
lg_a = md_a + ch*ceil((low_wl + low_wpml)/ch);

fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('omega/(2*pi) = %.2d,   1/h = %d   1/ch = %d,  NPW = %d \n',high_omega/(2*pi), 1/h, 1/ch, NPW);


%% Step 1: Solve the Hemholtz equation with the same source but with a relative low frequency sqrt(\omega) by Standard FEM, mesh size \omega*h = constant
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Step1: S-FEM, low frequency\n');

tic;
a = lg_a;                      % large computational domain [-lg_a, lg_a]^2
[node,elem] = squaremesh([-a,a,-a,a],h);
omega = low_omega;
wpml = low_wpml; sigmaMax = low_sigmaMax;
sigma = 1/100;
source = @(x) -4*sqrt(omega)*1i*1/(2*pi*sigma^2)*exp( -( (x(:,1)+1/4).^2 + (x(:,2)+1/4).^2  )/(2*sigma^2) );

A = assemble_Helmholtz_matrix_SFEM(node,elem,omega,wpml,sigmaMax,speed,fquadorder);
b = assemble_RHS(node,elem,source,fquadorder);
[~,bdEdge,isBdNode] = findboundary(elem);
freeNode = find(~isBdNode);
low_u = zeros(size(node(:,1)));
tic;
low_u(freeNode) = A(freeNode,freeNode)\b(freeNode);
toc;
FJ_showresult(node,elem,real(low_u));

%% Step 2: Use MPM to find ray directions d_c with low frequency sqrt(\omega)
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep2: MPM, low frequency\n');

a = md_a;                     % middle computational domain [-md_a, md_a]^2
[mnode,melem] = squaremesh([-a,a,-a,a],h);
[cnode,celem] = squaremesh([-a,a,-a,a],ch);
cN = size(cnode,1);
cnumray = zeros(cN,Nray);
cnumang = zeros(cN,Nray);
cray = ex_ray(cnode,xc,yc,1);

fprintf('MPM time: \n');
tic;
for i = 1:cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    d0 = sqrt((x0-xc)^2 + (y0-yc)^2);
    if d0 <= dr               % ray directions near source point
        cnumray(i,:) = cray(i,:);
    else                     % ray directions computed by NMLA far away from source point
        wl = 2*pi*speed([x0,y0])/low_omega;
        rh = wl/NPW;
        r = -wl:rh:wl;  r = r';
        
        est_ang0 = ex_ray([x0,y0],xc,yc,0);
        
        est_ang = est_ang0 + pi/6;
        x = r*cos(est_ang) + x0;  y = r*sin(est_ang) + y0;
        xy = [x, y];
        u = interpolation(node,elem,xy,low_u);
        u = u.*sqrt((x-xc).^2 + (y-yc).^2).^(1/2);
        [z] = Matrix_Pencil(u,1);
        
        dang1 = real(acos(1i*imag(log(z))/(1i*low_omega*rh)));
        ang1 = [est_ang + dang1, est_ang - dang1];
        ang1 = principal_angle(ang1);
        
        
        
        %% MPM sampling on one direction
        est_ang = est_ang0 + pi/3;
        x = r*cos(est_ang) + x0;  y = r*sin(est_ang) + y0;
        xy = [x, y];
        u = interpolation(node,elem,xy,low_u);
        u = u.*sqrt((x-xc).^2 + (y-yc).^2).^(1/2);
        [z] = Matrix_Pencil(u,1);
        
        dang2 = real(acos(1i*imag(log(z))/(1i*low_omega*rh)));
        ang2 = [est_ang + dang2, est_ang - dang2];
        ang2 = principal_angle(ang2);
        
        %% MPM sampling on another orthogonal direction
        est_ang =  est_ang0 + 3*pi/5;
        x = r*cos(est_ang) + x0;  y = r*sin(est_ang) + y0;
        xy = [x, y];
        u = interpolation(node,elem,xy,low_u);
        u = u.*sqrt((x-xc).^2 + (y-yc).^2).^(1/2);
        [z] = Matrix_Pencil(u,1);
        
        dang3 = real(acos(1i*imag(log(z))/(1i*low_omega*rh)));
        ang3 = [est_ang + dang3, est_ang - dang3];
        ang3 = principal_angle(ang3);
        
        %% correction
        temp = find_angle_3(ang1, ang2, ang3, eps);
        if (length(temp) == 1)
            cnumang(i) = temp;
        else
            tol = pi/360;
            temp = find_angle_3(ang1, ang2, ang3, tol);
            while (length(temp) ~= 1)  && (tol < 20*pi/180)
                tol = tol + pi/360;
                temp = find_angle_3(ang1, ang2, ang3, tol);
            end
            cnumang(i) = mean(temp);
        end
        cnumray(i) = exp(1i*cnumang(i));
        
    end
end
toc;

numray = interpolation(cnode,celem,mnode,cnumray);

ray = ex_ray(mnode,xc,yc,1);
md = sqrt((mnode(:,1)-xc).^2 + (mnode(:,2)-yc).^2);
ray = ray.*(1 - (md<eps));

% ray directions near source point
numray = numray.*(md>dr) + ray.*(md<=dr);
exray = ray;
dray = numray-exray;
norm(dray,inf)
figure(1);
FJ_showresult(mnode,melem,real(dray));


%% Step 3: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_c
% numray = exray;
node = mnode; elem = melem; a = md_a;
N = length(node(:,1));
omega = high_omega;
sigma = 1/100;
source = @(x) -4*sqrt(omega)*1i*1/(2*pi*sigma^2)*exp( -( (x(:,1)+1/4).^2 + (x(:,2)+1/4).^2  )/(2*sigma^2) );
[A,M] = assemble_Helmholtz_matrix_with_ray_1(node,elem,omega,wpml,sigmaMax,speed,numray,fquadorder);
b = assemble_RHS_with_ray_1(node,elem,omega,source,speed,ray,fquadorder);

% n = round(2*a/h + 1);
% xn = round((xc + a)/h);
% yn = round((yc + a)/h);
% xyn = n*xn + yn + 1;  
% b = M(:,xyn)/h^2;

[~,~,isBdNode] = findboundary(elem);
freeNode = find(~isBdNode);
v = zeros(N,1);
tic;
v(freeNode) = A(freeNode,freeNode)\b(freeNode);
toc;


%% reconstruct the solution
grad = numray(:);
grad = [real(grad),imag(grad)];
repnode = repmat(node,Nray,1);
temp = grad(:,1).*repnode(:,1) + grad(:,2).*repnode(:,2);

c = speed(node);    % medium speed
k = omega./c;           % wavenumber
kk = repmat(k,1,Nray);

u = v.*exp(1i*kk(:).*temp);
u = reshape(u,N,Nray);
u = sum(u,2);
u1 = u;
uh = u1;
% uh = u1/(sum(u1)*h^2);

figure(2);
FJ_showresult(node,elem,real(u1))


%% reference solution
NPW = 60;
rh = 1/2^(round(log2((NPW*high_omega)/(2*pi)))); 
% rh = 1/1600;

fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Reference S-FEM:\n omega/(2*pi) = %.2d,   1/rh = %d   NPW = %d \n',omega/(2*pi), 1/rh, NPW);
[rnode,relem] = squaremesh([-a,a,-a,a],rh);
[rA,rM] = assemble_Helmholtz_matrix_SFEM(rnode,relem,omega,wpml,sigmaMax,speed,fquadorder);
% rn = round(2*a/rh + 1);
% rxn = round((xc + a)/rh);
% ryn = round((yc + a)/rh);
% rxyn = rn*rxn + ryn + 1; 
% rb = rM(:,rxyn)/(1/2*rh^2);

rb = assemble_RHS(rnode,relem,source,fquadorder);

[~,~,isBdNode] = findboundary(relem);
freeNode = find(~isBdNode);
ru = zeros(size(rnode(:,1)));
ru(freeNode) = rA(freeNode,freeNode)\rb(freeNode);
% FJ_showresult(rnode,relem,real(ru));
%% Plotting
show_ray_solution(a,h,rh,uh,ru,0.45)





