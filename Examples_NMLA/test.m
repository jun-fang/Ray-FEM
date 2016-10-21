%% WT normalize the column of mass matrix
h = 1/128;
[node,elem] = squaremesh([-1,1,-1,1],h);
f = nodal_basis(0,0,1,node);
sum(f)*(h*h)



if (0)
add_omega = [120 140 160]*pi;
idx = [1 2 3 4 5 6 8 10];
rayerr2 = 0*add_omega; 
uerr2 = 0*add_omega;
rayerr1 = rayerr2;
uerr1 = uerr2;
for i = 1:length(add_omega)
    rayerr1(i) = test12(add_omega(i),-0.48,rec_omega,rec_ang_err1);
    uerr1(i) = test12(add_omega(i),-0.5,rec_omega,rec_NR_err1);
    rayerr2(i) = test12(add_omega(i),-0.48,rec_omega,rec_ang_err2);
    uerr2(i) = test12(add_omega(i),-0.5,rec_omega,rec_NR_err2);
end

omega = [rec_omega, add_omega];
rayerr1 = [rec_ang_err1, rayerr1];
uerr1 = [rec_NR_err1, uerr1];
rayerr2 = [rec_ang_err2, rayerr2];
uerr2 = [rec_NR_err2, uerr2];

omega = omega(idx);
rayerr1 = rayerr1(idx);
rayerr2 = rayerr2(idx);
uerr1 = uerr1(idx);
uerr2 = uerr2(idx);

save('result2_ex2.mat','omega','rayerr1','rayerr2','uerr1','uerr2');
end

if (0)
pde = Helmholtz_data1;
xs = 2; ys = 2;            % point source location

%% Set up
plt = 0;                   % show solution or not
fquadorder = 3;            % numerical quadrature order
solver = 'DIR';            % linear system solver
Nray = 1;                  % one ray direction
sec_opt = 1;               % NMLA second order correction or not

global omega a;
lg_a = 5/4;  md_a = 7/8;  sm_a = 1/2;
Rest = 1;  NPW = 6;

cp_omega = [40 80 160 320]*pi;
N = 1;
ang_err = zeros(1,N);
r = zeros(1,N);
for ii = 1:N
high_omega = cp_omega(ii);  low_omega = sqrt(high_omega);
h = 1/(NPW*round(high_omega/(2*pi)));
ch = 1/(2*NPW*round(low_omega/(2*pi)));

if (high_omega >= 100*pi) && (high_omega < 160*pi)
    lg_a = 1;  md_a = 2/3;  sm_a = 1/2;
elseif high_omega >= 160*pi
    lg_a = 7/8;  md_a = 5/8;  sm_a = 1/2;
end

a = md_a;  omega = high_omega;
[mnode,melem] = squaremesh([-a,a,-a,a],h);
uh1 = pde.ex_u(mnode);

%% Step 4: NMLA to find original ray directions d_o with wavenumber k
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nCase %d \n', ii);

a = sm_a;
[node,elem] = squaremesh([-a,a,-a,a],h);

[cnode,celem] = squaremesh([-a,a,-a,a],ch);
cN = size(cnode,1);
cnumray = zeros(cN,Nray);
[ux,uy] = num_derivative(uh1,h,2);
Du = pde.Du(mnode);
ux = Du(:,1);  uy = Du(:,2);

fprintf('NMLA time: \n');
tic;
for i = 1:cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    c0 = pde.speed(cnode(i,:));
    [cnumray(i,:),r(ii)] = NMLA_2D_correction(x0,y0,c0,omega,Rest,mnode,melem,uh1,ux,uy,pde,1/5,Nray,'num',sec_opt,plt);
end
numray2 = interpolation(cnode,celem,node,cnumray);
toc;

exray = ex_ray_angle(node,xs,ys);
diffang2 = numray2 - exray;
ang_err(ii) = h*norm(diffang2,2)/(h*norm(exray,2));
end
ang_err
% r

% showrate(cp_omega(1:N),ang_err)

end
