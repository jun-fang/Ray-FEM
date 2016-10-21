%% Reference solution for one point source problem: inside domain in homogeneous medium
clear;
xc = -3/8;   yc = -3/8;           % point source location
speed = @(x) ones(size(x(:,1)));
   
plt = 0;                         % plot solution or not
fquadorder = 3;                  % numerical quadrature order
Nray = 1;                        % no ray crosssing, number of ray direction is 1

omega = 60*pi;              % high frequency
NPW = 40;                         % number of grid points per wavelength

h = 1/2^(round(log2((8*omega)/(2*pi))));    %% 1/h should be a mutiple of 40!!
ch = 8*h;
rh = 1/2^(round(log2((NPW*omega)/(2*pi))));    %% 1/h should be a mutiple of 40!!
h = rh;

wl = 2*pi/omega;
wpml = 2*ch*ceil(wl/ch);     % width of PML
npml = round(ceil(wpml/h));
sigmaMax = 25/wpml;                  % Maximun absorbtion

fprintf('omega/(2*pi) = %d,   1/h = %d,   NPW = %d \n',omega/(2*pi), 1/rh, NPW);

a = 1/2 + ch*ceil(wpml/ch);
n = round(2*a/h + 1);
xn = round((xc + a)/h);
yn = round((yc + a)/h);
xyn = n*xn + yn + 1; 


%% Direct solver
if (0)
ttstart = cputime;
[node,elem] = squaremesh([-a,a,-a,a],h);
[A,M] = assemble_Helmholtz_matrix_SFEM(node,elem,omega,wpml,sigmaMax,speed,fquadorder);
b = M(:,xyn)/h^2;
[~,~,isBdNode] = findboundary(elem);
freeNode = find(~isBdNode);
ru = zeros(size(node(:,1)));
ru(freeNode) = A(freeNode,freeNode)\b(freeNode);
ttend = cputime;
fprintf('Direct solver time: %d\n', ttend -ttstart);

% norm(UGlobal - ru,inf)
% figure(2);
% FJ_showresult(node,elem,real(ru));
% clear node elem A b freeNode isBdNode;


u1 = ru; h1 = rh;
end





%% Direct solver
if (0)

NPW = 20;                         % number of grid points per wavelength

h = 1/2^(round(log2((8*omega)/(2*pi))));    %% 1/h should be a mutiple of 40!!
ch = 8*h;
rh = 1/2^(round(log2((NPW*omega)/(2*pi))));    %% 1/h should be a mutiple of 40!!
h = rh;

wl = 2*pi/omega;
wpml = 2*wl;     % width of PML
npml = round(ceil(wpml/h));
sigmaMax = 25/wpml;                  % Maximun absorbtion

fprintf('omega/(2*pi) = %d,   1/h = %d,   NPW = %d \n',omega/(2*pi), 1/rh, NPW);

a = 1/2 + ch*ceil(wpml/ch);
n = round(2*a/h + 1);
xn = round((xc + a)/h);
yn = round((yc + a)/h);
xyn = n*xn + yn + 1; 



ttstart = cputime;
[node,elem] = squaremesh([-a,a,-a,a],h);
[A,M] = assemble_Helmholtz_matrix_SFEM(node,elem,omega,wpml,sigmaMax,speed,fquadorder);
b = M(:,xyn)/h^2;
[~,~,isBdNode] = findboundary(elem);
freeNode = find(~isBdNode);
ru = zeros(size(node(:,1)));
ru(freeNode) = A(freeNode,freeNode)\b(freeNode);
ttend = cputime;
fprintf('Direct solver time: %d\n', ttend -ttstart);

% norm(UGlobal - ru,inf)
% figure(2);
FJ_showresult(node,elem,real(ru));
clear node elem A b freeNode isBdNode;

u2 = ru; h2 = rh;


%% plot
figure(3);
n1 = round(2*a/h1) + 1;
uu1 = reshape(u1,n1,n1);

n2 = round(2*a/h2) + 1;
uu2 = reshape(u2,n2,n2);

hold off;
dy = 3/4;
yn = round(dy/h) + 1;
xx = -a:h:a;
uh2 = uu2(yn,:);

yn1 = round(dy/h1) + 1;
xx1 = -a:h1:a;
uh1 = uu1(yn1,:);
plot(xx,real(uh2),'ro-');
hold on;
plot(xx1,real(uh1),'k');
xlabel('x');
ylabel('Wavefield');
legend('Ray-FEM solution','Reference solution','LOCATION','Best');



end









%% Leo's solver
if (1)
ttstart=cputime;

[node,elem] = squaremesh([-a,a,-a,a],h);

% Number of subdomains
nSub = 5;
% initialazing the  global model 
M0 = model;
M0.init(node,elem,omega,wpml,sigmaMax,speed,fquadorder);

% % performing LU factorization (only if we want ot compare against the true
% % solution 
% M0.LU_factorization()

% solving the PDE using the global model 
% f = assemble_RHS(node,elem,source,fquadorder);
[~,M] = assemble_Helmholtz_matrix_SFEM(node,elem,omega,wpml,sigmaMax,speed,fquadorder);

f = M(:,xyn)/h^2;

% getting only the interior degrees of freedom
fInt = f(M0.freeNode);

%% Local problems ( We need to find a way to perform the partitioning)
% This needs to be properly encapsulated
% geometrical information 
x = unique(node(:,1));
y = unique(node(:,2));

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
                                                min(y),max(y)],h);
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
                    wpmlvec,sigmaMax,speed,fquadorder);

end

%% we need to define all the local (and global) local Indices
% This should be encapsulated; however, given that the meshes can be
% arbitrary we need to do it by had.

% We start by defining the indices associated to each boundary (locally)
indxn = {};
indxnp = {};
indx0  = {};
indx1  = {};

for ii = 1:nSub
    if ii ~= 1
        indx0{ii} = find(MArray{ii}.node(MArray{ii}.freeNode,1) == separatorN(ii-1) );
        indx1{ii} = find(MArray{ii}.node(MArray{ii}.freeNode,1) == separator1(ii-1) );
    end
    
    if ii ~= nSub
        indxn{ii} = find(MArray{ii}.node(MArray{ii}.freeNode,1) == separatorN(ii) );
        indxnp{ii}= find(MArray{ii}.node(MArray{ii}.freeNode,1) == separator1(ii) );
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
        indIntGlobal{ii} = find(M0.node(M0.freeNode,1) >= separator1(ii-1) );
    else
        indIntGlobal{ii} = find((M0.node(M0.freeNode,1) <= separatorN(ii) ).*  ...
                                (M0.node(M0.freeNode,1) >= separator1(ii-1)));
    end
end
  
% building the local-local indices
indIntLocal = {};
for ii = 1:nSub
    if ii == 1
        indIntLocal{ii} = find(MArray{ii}.node(MArray{ii}.freeNode,1) <= separatorN(ii) );
    elseif ii == nSub
        indIntLocal{ii} = find(MArray{ii}.node(MArray{ii}.freeNode,1) >= separator1(ii-1) );
    else
        indIntLocal{ii} = find((MArray{ii}.node(MArray{ii}.freeNode,1) <= separatorN(ii) ).*  ...
                                (MArray{ii}.node(MArray{ii}.freeNode,1) >= separator1(ii-1)));
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
 
UGlobal = zeros(size(f,1),1);
UGlobal(M0.freeNode) = U;

ttend=cputime;

fprintf('Leo''s solver time: %d\n', ttend - ttstart);

clear node elem UsolArray elemArray fInt indIntGlobal indIntLocal;
clear nodeArray f fGammaPol MMu U fGamma uArray uGmres; 
ru = UGlobal;
% figure(1); 
% M0.showresult(real(UGlobal));

% norm(u1 - ru,inf)
end

save('ex1_reference_solution.mat','rh','ru','a','wpml','sigmaMax');

