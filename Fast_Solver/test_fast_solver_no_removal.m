%% One point source inside domain (homogeneous medium):  iterative idea
% script that computes the low and high frequency systems using the 
% method of polarized traces

addpath(genpath('../../ifem/'));
addpath('../Methods/');
addpath('../NMLA/');
addpath('../Functions/')
addpath('../Plots_Prints/');
addpath(genpath('../Classes/'));
addpath('../Singularity_Treatment/');


xs = 0; ys = 0;                     % source location 
epsilon = 1/(2*pi);                 % cut-off parameter   speed = @(x) ones(size(x,1),1);    % medium speed
speed = @(x) ones(size(x,1),1);      % medium speed



plt = 1;                           % plot solution or not
fquadorder = 6;                    % numerical quadrature order
Nray = 1;
Rest = 1;
pde = [];
pct = 1/5;
data = 'num';
opt = 1;

high_omega = 60*pi;
low_omega = sqrt(high_omega);
NPW = 10;

% number of subdomains
nSub = round(high_omega / (10*pi));

h = 1/(20*round(high_omega*NPW/(2*pi*20)));
% high frequency h
hTilde = 1/(20*max(round(low_omega*NPW/(2*pi*20)),1));
% low frequency h (tildeh)

wavelength = 2*pi/high_omega;  
wpml = ceil(wavelength/h)*h;               % width of PML
sigmaMax = 25/wpml;                % Maximun absorbtion

% observation radious
r = 8*wpml;

lg_a = 1;
md_a = 0.65;

sm_a = 1/2;

fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('omega/(2*pi) = %d,   1/h = %d   1/ch = %d,  NPW = %d \n',high_omega/(2*pi), 1/h, 1/hTilde, NPW);

sigma = 1/100;
k = high_omega;


%% Step 1: Solve the Hemholtz equation with the same source but with a relative low frequency sqrt(\omega) by Standard FEM, mesh size \omega*h = constant
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('Step1: S-FEM, low frequency\n');

a = lg_a;
[lnode,lelem] = squaremesh([-a,a,-a,a],h);

fprintf('Starting the construction of the fast solver \n');
omega = low_omega;

M0 = model;
M0.init(lnode,lelem,low_omega,wpml,h/NPW, speed, fquadorder);

npml = wpml/h;

fquadorder = 4;  
% assembling the RHS 
%f = assemble_RHS(nodeLow,elemLow,source,fquadorder);
f = assemble_RHS_SFEM_with_ST(lnode,lelem,xs,ys,omega,wpml,sigmaMax,epsilon,fquadorder);


% this should only be for comparison purpouses 
% factorizing the matrix within the model class
% M0.LU_factorization()
% solving the PDE using the global model 
% u_std = M0.solve(f);
    


% getting only the interior degrees of freedom
fInt = f(M0.freeNode);

%[u_std] = Standard_FEM_PML_PointSource(lnode,lelem,omega,wpml,sigmaMax,xs,ys,speed,fquadorder,plt);



%% Local problems ( We need to find a way to perform the partitioning)
% This needs to be properly encapsulated
% geometrical information 
% x = unique(nodeLow(:,1));
% y = unique(nodeLow(:,2));

x = unique(M0.node(:,1));
y = unique(M0.node(:,2));


% gotta be careful with the PML (they have to be outside the 
% interior domain
xInt = x(npml:end-npml);
nInt = length(xInt);

% number of points per subdomains
nSubInt = round(nInt/nSub);

% the limits for the domaind decomposition
nIndLim = round(linspace(1,nInt-1,nSub+1));

% indices at which each subdomain starts (and finishes)
indn = npml + nIndLim(2:end);
ind1 = npml + nIndLim(1:end-1)+1;
ind1(1) = ind1(1) -1;

% last and first coordinate
xn = x(indn);
x1 = x(ind1);

nodeArray = {};
elemArray = {};
for ii = 1:nSub
    [nodeArray{ii},elemArray{ii}] = squaremesh([x(ind1(ii)-npml),...
                                                x(indn(ii)+npml),...
                                                -a,a],h);
end

% we define the physical location of the interface boundaries
separatorN = xn(1:end-1);
separator1 = x1(2:end);

% building the array of objects
MArray = {};
SubArray = {};
for ii = 1:nSub
    % we will need to make sure that the wpml has the correct width
    MArray{ii} = model;
    % we need to be carefull when defining the PML
    if ii == 1
        wpmlvec = [wpml, wpml-2*h, wpml, wpml];
    elseif ii == nSub
        wpmlvec = [wpml-2*h, wpml, wpml, wpml];
    else
        wpmlvec = [wpml-2*h, wpml-2*h, wpml, wpml];
    end     
    MArray{ii}.init(nodeArray{ii},elemArray{ii},omega,...
                    wpmlvec,h/NPW,speed,fquadorder);

end

%% we need to define all the local (and global) local Indices
% We start by defining the indices associated to each boundary (locally)
indxn = {};
indxnp = {};
indx0  = {};
indx1  = {};

for ii = 1:nSub
    if ii ~= 1
        indx0{ii} = find( abs(MArray{ii}.node(MArray{ii}.freeNode,1) - separatorN(ii-1)) < h/10 );
        indx1{ii} = find( abs(MArray{ii}.node(MArray{ii}.freeNode,1) - separator1(ii-1)) < h/10 );
    end
    
    if ii ~= nSub
        indxn{ii} = find( abs(MArray{ii}.node(MArray{ii}.freeNode,1) - separatorN(ii)) < h/10 );
        indxnp{ii}= find( abs(MArray{ii}.node(MArray{ii}.freeNode,1) - separator1(ii)) < h/10 );
    end
end
indxn{nSub} = [];
indxnp{nSub} = [];

% building the the local-global indices

indIntGlobal = {};
for ii = 1:nSub
    if ii == 1
        indIntGlobal{ii} = find(M0.node(M0.freeNode,1) <= separatorN(ii) );
    elseif ii == nSub
        indIntGlobal{ii} = find(M0.node(M0.freeNode,1) > separator1(ii-1)-h/10 );
    else
        indIntGlobal{ii} = find((M0.node(M0.freeNode,1) <= separatorN(ii) + h/10 ).*  ...
                                (M0.node(M0.freeNode,1) > separator1(ii-1) - h/10 ));
    end
end
  
 
% building the local-local indices

indIntLocal = {};
for ii = 1:nSub
    if ii == 1
        indIntLocal{ii} = find(MArray{ii}.node(MArray{ii}.freeNode,1) <= separatorN(ii) );
    elseif ii == nSub
        indIntLocal{ii} = find(MArray{ii}.node(MArray{ii}.freeNode,1) > separator1(ii-1)-h/10  );
    else
        indIntLocal{ii} = find((MArray{ii}.node(MArray{ii}.freeNode,1) <= separatorN(ii)+h/10 ).*  ...
                                (MArray{ii}.node(MArray{ii}.freeNode,1) > separator1(ii-1)-h/10 ));
    end
end


for ii = 1:nSub
    if size(indIntLocal{ii},1) ~= size(indIntGlobal{ii},1)
        fprintf('Mismatch of number of indices ot %d th layer \n', ii);
    end
end



%% building the subdomains from the models and the local indices
for ii = 1:nSub
    if ii == 1
        position = 'N';
    elseif ii == nSub
        position = 'S';
    else
        position = 'C'; %center
    end
    
    SubArray{ii} = subdomain;    
    SubArray{ii}.init_Model(MArray{ii},indIntLocal{ii}, indIntGlobal{ii}, indxn{ii}, indxnp{ii}, indx0{ii}, indx1{ii}, position)
end

%% Factorizing the local matrices
tic();
for ii = 1:nSub
   SubArray{ii}.LU_factorization();
end
toc();

%% bulding the global array 
DDM = layer_stripping;

DDM.init_model(SubArray,M0, nSub);

%% Solving local problems

uArray =  DDM.solve_local(fInt);

fGamma = DDM.formRHS(uArray);

%% Testing the application of the matrix M
% 
% DDM.MBase.LU_factorization();
% uGlobal = DDM.MBase.solveInt(fInt);
% uGamma = DDM.extractGlobalTraces(uGlobal);
%  
% res =  DDM.applyM(uGamma) + fGamma;
% 
% norm(res)

%% testing the application of the  matrix MM
fGammaPol = DDM.formPolarizedRHS(uArray);

MMu = DDM.applyMpolarized(fGammaPol );

DDM.initPermutationMatrix();
MM = @(u) DDM.applyMpolarized(u);

%% testing the Dinv 
LL = @(u) DDM.applyGSinv(u);

tic();
uGmres = gmres(MM, -fGammaPol,[], 1e-9, 300, LL );
toc()
%% testing that the solution coincides

u = uGmres(1:end/2) + uGmres(end/2+1:end);

res = DDM.applyMpolarized(uGmres) + fGammaPol;
norm(res)


%% reconstructing the solution within the subdomain
UsolArray = DDM.solveLocal(fInt, u);
U = DDM.reconstruct(UsolArray);
 
u_std = zeros(size(f,1),1);
u_std(M0.freeNode) = U;

% plotting the global smooth solution 
if plt == 1
    M0.showresult(real(u_std));
end

%% reconstruction of the full wave field

    % singular part
    x = lnode(:,1); y = lnode(:,2);
    rr = sqrt((x-xs).^2 + (y-ys).^2);
    ub = 1i/4*besselh(0,1,omega*rr);
    cf = cutoff(epsilon,2*epsilon,lnode,xs,ys);
    
    % low frequency solution: smooth + singularity
    u_low = u_std + ub.*cf;
    toc;
    
    
% plotting the global solution 
if plt == 1
    M0.showresult(real(u_low));
end  


%% Step 2: Use NMLA to find ray directions d_c with low frequency sqrt(\omega)
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep2: NMLA, low frequency\n');

[ux,uy] = num_derivative(u_std,h,2);

a = md_a;
[mnode,melem] = squaremesh([-a,a,-a,a],h);
[cnode,celem] = squaremesh([-a,a,-a,a],hTilde);
cN = size(cnode,1);
cnumray = zeros(cN,Nray);
cray = ex_ray(cnode,xs,ys,0);
cr1 = zeros(cN,1);
cd1 = cr1;

fprintf('NMLA time: \n');
tic;
for i = 1:cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    if d0 <= r
        cnumray(i,:) = cray(i,:);
    else
        Rest = d0;
        if d0 <= 2*r
            Rest = 2*d0;
        end
        c0 = speed(cnode(i,:));
        %[cnumray(i,:),~,cr1(i)] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,nodeLow,elemLow,u_std,ux,uy,pde,pct,Nray,data,opt,plt);
        [cnumray(i,:), ~] = NMLA(x0,y0,c0,omega,Rest,lnode,lelem,u_std,ux,uy,pde,pct,Nray,data,opt,plt);
        cd1(i) = cr1(i) - d0;
%         if r1 > d0
%             cr1(i) = 1;
%         end
    end
end
toc;

clear lnode lelem;

cdiffang = angle_error(cnumray,cray);
norm(cdiffang,2)/norm(cray);
norm(cdiffang,inf)

cnumray = exp(1i*cnumray);
numray1 = interpolation(cnode,celem,mnode,cnumray);

ray = ex_ray(mnode,xs,ys,0);
ray = exp(1i*ray);
md = sqrt((mnode(:,1)-xs).^2 + (mnode(:,2)-ys).^2);
ray = ray.*(1 - (md<eps));

numray1 = numray1.*(md>r) + ray.*(md<=r);
diffray1 = numray1 - ray;


%% Step 3: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_c
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep3: Ray-FEM, high frequency \n');

omega = high_omega;

MHigh = model;
MHigh.initRay(mnode,melem,high_omega,wpml,h/NPW, speed, numray1, fquadorder);

% solving the PDE using the global model 

option = 'homogeneous';
fHigh = assemble_RHS_RayFEM_with_ST(mnode,melem,xs,ys,omega,epsilon,...
                                    wpml,sigmaMax,ray,speed,fquadorder,option);

% % performing LU factorization (only if we want ot compare against the true
% % solution 
% MHigh.LU_factorization()
% % solving the PDE using the Multifrontal method
% u1 = MHigh.solve(fHigh);

fInt = fHigh(MHigh.freeNode);


%% starting to build the fast solver for the high frequency solver

fprintf('Building the data structure for the fast solver \n')

x = unique(MHigh.node(:,1));
y = unique(MHigh.node(:,2));

% gotta be careful with the PML (they have to be outside the 
% interior domain
xInt = x(npml:end-npml);
nInt = length(xInt);

% number of points per subdomains
nSubInt = round(nInt/nSub);

% the limits for the domaind decomposition
nIndLim = round(linspace(1,nInt-1,nSub+1));

% indices at which each subdomain starts (and finishes)
indn = npml + nIndLim(2:end);
ind1 = npml + nIndLim(1:end-1)+1;
ind1(1) = ind1(1) -1;

% we need to add this to the ray construction

% last and first coordinate
xn = x(indn);
x1 = x(ind1);

nodeArray = {};
elemArray = {};
for ii = 1:nSub
    [nodeArray{ii},elemArray{ii}] = squaremesh([x(ind1(ii)-npml),...
                                                x(indn(ii)+npml),...
                                                -a,a],h);
end

% we define the physical location of the interface boundaries
separatorN = xn(1:end-1);
separator1 = x1(2:end);

nx = size(x,1);
rayArray = {};
rayMatrix = (reshape(numray1,nx,nx)).';

for ii = 1:nSub
     aux = (rayMatrix(ind1(ii)-npml:indn(ii)+npml,:)).';
     rayArray{ii} = aux(:);
end


% building the array of objects
MArray = {};
SubArray = {};
for ii = 1:nSub
    % we will need to make sure that the wpml has the correct width
    MArray{ii} = model;
    % we need to be carefull when defining the PML
    if ii == 1
        wpmlvec = [wpml, wpml-2*h, wpml, wpml];
    elseif ii == nSub
        wpmlvec = [wpml-2*h, wpml, wpml, wpml];
    else
        wpmlvec = [wpml-2*h, wpml-2*h, wpml, wpml];
    end     
    MArray{ii}.initRay(nodeArray{ii},elemArray{ii},high_omega,...
                    wpmlvec,h/NPW,speed,rayArray{ii},fquadorder);

end

%% we need to define all the local (and global) local Indices
% We start by defining the indices associated to each boundary (locally)
indxn = {};
indxnp = {};
indx0  = {};
indx1  = {};

for ii = 1:nSub
    if ii ~= 1
        indx0{ii} = find( abs(MArray{ii}.node(MArray{ii}.freeNode,1) - separatorN(ii-1)) < h/10 );
        indx1{ii} = find( abs(MArray{ii}.node(MArray{ii}.freeNode,1) - separator1(ii-1)) < h/10 );
    end
    
    if ii ~= nSub
        indxn{ii} = find( abs(MArray{ii}.node(MArray{ii}.freeNode,1) - separatorN(ii)) < h/10 );
        indxnp{ii}= find( abs(MArray{ii}.node(MArray{ii}.freeNode,1) - separator1(ii)) < h/10 );
    end
end
indxn{nSub} = [];
indxnp{nSub} = [];

% building the the local-global indices

indIntGlobal = {};
for ii = 1:nSub
    if ii == 1
        indIntGlobal{ii} = find(MHigh.node(MHigh.freeNode,1) <= separatorN(ii) + h/10  );
    elseif ii == nSub
        indIntGlobal{ii} = find(MHigh.node(MHigh.freeNode,1) > separator1(ii-1) - h/10  );
    else
        indIntGlobal{ii} = find((MHigh.node(MHigh.freeNode,1) <= separatorN(ii) + h/10 ).*  ...
                                (MHigh.node(MHigh.freeNode,1) > separator1(ii-1)- h/10 ));
    end
end
  
 
% building the local-local indices

indIntLocal = {};
for ii = 1:nSub
    if ii == 1
        indIntLocal{ii} = find(MArray{ii}.node(MArray{ii}.freeNode,1)  <= separatorN(ii) + h/10  );
    elseif ii == nSub
        indIntLocal{ii} = find(MArray{ii}.node(MArray{ii}.freeNode,1) > separator1(ii-1) - h/10 );
    else
        indIntLocal{ii} = find((MArray{ii}.node(MArray{ii}.freeNode,1) <= separatorN(ii)+h/10 ).*  ...
                                (MArray{ii}.node(MArray{ii}.freeNode,1) > separator1(ii-1)-h/10 ));
    end
end


for ii = 1:nSub
    if size(indIntLocal{ii},1) ~= size(indIntGlobal{ii},1)
        fprintf('Mismatch of number of indices ot %d th layer \n', ii);
    end
end



%% building the subdomains from the models and the local indices
for ii = 1:nSub
    if ii == 1
        position = 'N';
    elseif ii == nSub
        position = 'S';
    else
        position = 'C'; %center
    end
    
    SubArray{ii} = subdomain;    
    SubArray{ii}.init_Model(MArray{ii},indIntLocal{ii}, indIntGlobal{ii}, indxn{ii}, indxnp{ii}, indx0{ii}, indx1{ii}, position)
end

%% Factorizing the local matrices
tic();
for ii = 1:nSub
   SubArray{ii}.LU_factorization();
end
toc();

%% bulding the global array 
DDMHigh = layer_stripping;

DDMHigh.init_model(SubArray,MHigh, nSub);

%% Solving local problems

uArray =  DDMHigh.solve_local(fInt);

fGamma = DDMHigh.formRHS(uArray);

%% Testing the application of the matrix M
% 
% DDM.MBase.LU_factorization();
% uGlobal = DDM.MBase.solveInt(fInt);
% uGamma = DDM.extractGlobalTraces(uGlobal);
%  
% res =  DDM.applyM(uGamma) + fGamma;
% 
% norm(res)

%% testing the application of the  matrix MM
fGammaPol = DDMHigh.formPolarizedRHS(uArray);

MMu = DDMHigh.applyMpolarized(fGammaPol );

DDMHigh.initPermutationMatrix();
MM = @(u) DDMHigh.applyMpolarized(u);

%% testing the Dinv 
LL = @(u) DDMHigh.applyGSinv(u);

tic();
uGmres = gmres(MM, -fGammaPol,[], 1e-6, 30, LL );
toc()
%% testing that the solution coincides

uBdy = uGmres(1:end/2) + uGmres(end/2+1:end);

res = DDMHigh.applyMpolarized(uGmres) + fGammaPol;
norm(res)


%% reconstructing the solution within the subdomain
UsolArray = DDMHigh.solveLocal(fInt, uBdy);
U = DDMHigh.reconstruct(UsolArray);
 
u1 = zeros(size(fHigh,1),1);
u1(MHigh.freeNode) = U;


N = size(mnode,1);       % number of grid points
Nray = size(numray1,2);     % number of rays crossing at each grid node
Ndof = N*Nray; 


if  plt == 1
    
    grad = numray1(:);
    grad = [real(grad),imag(grad)];
    repnode = repmat(mnode,Nray,1);
    temp = grad(:,1).*repnode(:,1) + grad(:,2).*repnode(:,2);

    c = speed(mnode);    % medium speed
    k = omega./c;       % wavenumber
    kk = repmat(k,1,Nray);

    u = u1.*exp(1i*kk(:).*temp);
    u = reshape(u,N,Nray);
    u = sum(u,2);

    % showing the real part of the solution
    figure(2); 
    MHigh.showresult(real(u));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% until here is should be OK


%% Step 4: NMLA to find original ray directions d_o with wavenumber k
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep4: NMLA, high frequency \n');

a = sm_a;
[node,elem] = squaremesh([-a,a,-a,a],h);
[cnode,celem] = squaremesh([-a,a,-a,a],hTilde);
cN = size(cnode,1);
cnumray = zeros(cN,Nray);

cray = ex_ray(cnode,xs,ys,0);
cr2 = zeros(cN,1);
cd2 = cr2;

[ux,uy] = num_derivative(u1,h,2);

fprintf('NMLA time: \n');
tic;
for i = 1:cN
    x0 = cnode(i,1);  y0 = cnode(i,2);
    d0 = sqrt((x0-xs)^2 + (y0-ys)^2);
    if d0 <= r
        cnumray(i,:) = cray(i,:);
    else
        Rest = d0;
        if d0 <= 2*r
            Rest = 2*d0;
        end
        c0 = speed(cnode(i,:));
        [cnumray(i,:)] = NMLA_2D_2nd(x0,y0,c0,omega,Rest,mnode,melem,u1,ux,uy,pde,pct,Nray,data,opt,plt);
        cd2(i) = cr2(i) - d0;
%         if cr > d0
%             cr2(i) = 1;
%         end
    end
end
toc;

% clear mnode melem;

cdiffang = angle_error(cnumray,cray);
norm(cdiffang,2)/norm(cray)
norm(cdiffang,inf)

cnumray = exp(1i*cnumray);
numray2 = interpolation(cnode,celem,node,cnumray);

ray = ex_ray(node,xs,ys,0);
ray = exp(1i*ray);
d = sqrt((node(:,1)-xs).^2 + (node(:,2)-ys).^2);
ray = ray.*(1 - (d < eps) );

numray2 = numray2.*(d>r) + ray.*(d<=r);
diffray2 = numray2 - ray;
% figure(2);
% FJ_showresult(node,elem,real(diffray2));

numray_dir = [real(numray2), imag(numray2)];
numray_dir = atan2(numray_dir(:,2), numray_dir(:,1));
numray_dir = numray_dir + 2*pi*(numray_dir<0);

diffang = angle_error(numray_dir,ex_ray(node,xs,ys,0));
diffang = diffang.*(d>r);
% figure(3);
% FJ_showresult(node,elem,real(diffang));
% title('NMLA angle error');
figure(3);
FJ_showresult(cnode,celem,cr2);

%% Step 5: Solve the original Helmholtz equation by Ray-based FEM with ray directions d_o
fprintf(['\n' '-'*ones(1,80) '\n']);
fprintf('\nStep5: Ray-FEM, high frequency \n');
% 
% omega = high_omega;
% [u2,~,v2] = Ray_FEM_PML_1_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,numray2,fquadorder,plt);
% 
% 
% %% 
% % [X,Y] = meshgrid(-a:h:a,-a:h:a);
% % [m,n] = size(X);
% % uh = reshape(u2,m,n);
% % 
% % figure(4)
% % contour(X,Y,real(uh));
% % title('Level set of Ray-FEM solution u_{ray}');
% 
% 
% %% map to polar
% figure(5);
% % r1 = 10.3*wpml;
% % r2 = 12.5*wpml;
% r1 = 0.834;
% r2 = 0.834 + 2.1*wavelength;
% theta1 = pi/4 - pi/16;
% theta2 = pi/4 + pi/16;
% subplot(2,1,1);
% mapto_polar(node,elem,omega,speed,v2,numray2,xs,ys,r1,r2,theta1,theta2);
% 
% subplot(2,1,2);
% mapto_polar(node,elem,omega,speed,v2,numray2,xs,ys,r1,r2,theta1,theta2);
% axis equal
% 
% 
% 
% %% 
% % figure(6);
% % % r1 = 0.81;
% % % r2 = 0.81 + 2.1*wavelength;
% % % omega = 80*pi;
% % [X,Y] = meshgrid(-a:1/1000:a,-a:1/1000:a);
% % [m,n] = size(X);
% % xynode = [X(:),Y(:)];
% % uu = ray_solution(node,elem,omega,speed,v2,numray2,xynode);
% % uu = reshape(uu,m,n);
% % mesh(X,Y,real(uu));
% % az = 0;
% % el = 90;
% % view(az, el);
% % axis equal; axis tight;
% % hold on;
% % 
% % theta = theta1:1/10000:theta2;
% % rr = r1: 1/5000:r2;
% % x1 = r1*cos(theta) + xs;   y1 = r1*sin(theta) + ys;
% % x2 = r2*cos(theta) + xs;   y2 = r2*sin(theta) + ys;
% % x3 = rr*cos(theta1) + xs;  y3 = rr*sin(theta1) + ys;
% % x4 = rr*cos(theta2) + xs;  y4 = rr*sin(theta2) + ys;
% % p = plot(x1,y1,'r-');
% % p.LineWidth = 3;  hold on;
% % p = plot(x2,y2,'r-');
% % p.LineWidth = 3;  hold on;
% % p = plot(x3,y3,'r-');
% % p.LineWidth = 3;  hold on;
% % p = plot(x4,y4,'r-');
% % p.LineWidth = 3;  hold on;
% % xlabel('x');
% % ylabel('y');

