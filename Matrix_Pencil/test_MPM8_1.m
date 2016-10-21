%% One point source problem:  iterative idea
% MPM
% outside domain (homogeneous medium): 
global omega;
global a;
global xs ys;

pde = Helmholtz_data1_0;

xs = 10;   ys = 10;             % point source location

plt = 0;                           % plot solution or not
fquadorder = 3;                    % numerical quadrature order
solver = 'DIR';
Nray = 1;





high_omega = 400*pi;
low_omega = sqrt(high_omega); %round(sqrt(high_omega)/4)*2*pi;
NPW = 8;


% h = 1/400; ch = 1/50;
h = 1/round((NPW*high_omega)/(2*pi));
ch = 1/(20*max(round(low_omega*NPW/(2*pi*20)),1));

high_wl = 2*pi/high_omega;  
low_wl = 2*pi/low_omega;  

sm_a = 1/2;
md_a = sm_a + high_wl;
lg_a = md_a + low_wl;
md_a = sm_a;
lg_a = sm_a + low_wl;

fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('omega/(2*pi) = %d,   1/h = %d   1/ch = %d,  NPW = %d \n',high_omega/(2*pi), 1/h, 1/ch, NPW);


%% Step 1: Solve the Hemholtz equation with the same source but with a relative low frequency sqrt(\omega) by Standard FEM, mesh size \omega*h = constant
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Step1: S-FEM, low frequency\n');

a = lg_a;
[lnode,lelem] = squaremesh([-a,a,-a,a],h);
omega = low_omega;
% [u_std,~,~,~] = Standard_FEM_IBC(lnode,lelem,omega,pde,fquadorder,solver,plt);

u_std = sqrt(omega)*besselh(0,1,omega*sqrt((lnode(:,1)-xs).^2 + (lnode(:,2)-ys).^2));
% FJ_showresult(lnode,lelem,real(u_std));

%% Step 2: Use MPM to find ray directions d_c with low frequency sqrt(\omega)
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep2: MPM, low frequency\n');


a = md_a;
[mnode,melem] = squaremesh([-a,a,-a,a],h);
[cnode,celem] = squaremesh([-a,a,-a,a],ch);
cN = size(cnode,1);
cnumray = zeros(cN,Nray);
cray = ex_ray(cnode,xs,ys,0);
exray = cray;
err1 = zeros(cN,Nray);
err2 = err1;


wl = low_wl;
rh = wl/NPW;

fprintf('MPM time: \n');
tic;
for i = 1:cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    r = -wl:rh:wl;
            r = r';
            %             est_ang =  exray(i)/1.1 + pi/12;
            est_ang =  exray(i);% + (-1)^(exray(i)>pi)*pi/30;
            x = r*cos(est_ang) + x0;
            y = r*sin(est_ang) + y0;
            xy = [x, y];
            
            nu = interpolation(lnode,lelem,xy,u_std);
            u = nu;
            u = u.*sqrt((x-xs).^2 + (y-ys).^2).^(1/2);
            
            %% MPM
            [z] = Matrix_Pencil_2(u);

            err1(i,:) = min(abs(abs(exray(i) - est_ang) - abs(acos(1i*imag(log(z))/(1i*omega*rh))) ));
            err2(i,:) = min(abs(cos(exray(i) - est_ang) - 1i*imag(log(z))/(1i*omega*rh)));
            
end
toc;
norm(err1,inf)
norm(err1)/norm(cray)

