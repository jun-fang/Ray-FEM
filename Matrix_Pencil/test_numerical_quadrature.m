%% Test for numerical quadrature

global omega xs ys a;
xs = 10; ys = 10;                 % source location
aa = 1/2;                         % computational domain [-1/2,1/2]^2

pde = MPM_data1;                  % one point source pde structure
Nray = 4;                         % one ray direction at every grid point
plt = 0;                          % plot solution or not
fquadorder = 3;                   % numerical quadrature order
solver = 'DIR';                   % direct method to solve Au=b by u=A\b

%% record frequencies and errors
rec_omega = [10 20 40 80 160]*pi;
rec_err1 = zeros(length(rec_omega),1);
rec_h = rec_err1;

nn = 1;

for j = 1:nn%length(rec_omega)
    tic;
    j
    high_omega = rec_omega(j);         % high frequency
    low_omega = sqrt(high_omega);      % low frequency
    omega = high_omega;
    NPW = 6;                           % number of grid points per wavelength
    h = 1/round((NPW*omega)/(2*pi));   % mesh size
    %     h = h/2^(j-2);
    %     rec_h(j) = h;
    %     h = 1/100;
    
    
    a = 1/2;
    [node,elem] = squaremesh([-a,a,-a,a],h);
    exray = pde.ray(node);
    
    % Exact Ray-FEM
    ray1 = exp(1i*exray);
    %     ray1 = exp(1i*rand(size(exray)));
    ray2 = exray - pi/4;
    ray2 = exp(1i*ray2);
%         ray = [ray1, ray2];
    ray = ray1;
    ray(10) = ray2(10);
    %     [uh,A,~,~,rel_L2_err] = Ray_FEM_IBC_1(node,elem,omega,pde,ray,fquadorder,plt);
    [~,A3,~,~,~] = Ray_FEM_IBC_1(node,elem,omega,pde,ray,3,plt);
    [~,A9,~,~,~] = Ray_FEM_IBC_1(node,elem,omega,pde,ray,9,plt);
    d = A3-A9;
    norm(d(:),inf)
    
    N = size(node,1);       % number of grid points
    NT = size(elem,1);      % number of triangle elements
    Nray = size(ray,2);     % number of rays crossing at each grid node
    Ndof = N*Nray;          % degree of freedom
    
    c = pde.speed(node);    % medium speed
    k = omega./c;           % wavenumber
    kk = repmat(k,1,Nray);
    
    
    %% %% 3rd order Numerical Quadrature
    fquadorder = 3;     
    [lambda,weight] = quadpts(fquadorder);
    phi = lambda;           % the value of piesewise linear basis functions at numerical quadrature nodes
    nQuad = size(lambda,1);
    
    
    %% Compute geometric quantities and gradient of local basis
    [Dphi,area] = gradbasis(node,elem);
    reparea = repmat(area,Nray,1);
    
    
    %% Assemble the linear system Av = b
    
    
    
    A = sparse(Ndof,Ndof);
    A1 = A;
    A2 = A;
    A3 = A;
    M = A;
    
    bt = zeros(NT*Nray,3);
    b = zeros(Ndof,1);
    
    % fprintf('Assembling time:\n');
    % tic;
    for p = 1:nQuad
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        reppxy = repmat(pxy,Nray,1);
        fp = pde.f(reppxy);
        k2 = (omega./pde.speed(pxy)).^2;
        
        for i = 1:3
            gradtempi = - ray(elem(:,i),:);
            gradtempi = gradtempi(:);
            gradtempi = [real(gradtempi), imag(gradtempi)];
            fphasei = gradtempi(:,1).*reppxy(:,1) + gradtempi(:,2).*reppxy(:,2);
            kki = kk(elem(:,i),:);
            kki = kki(:);
            phasei = exp(1i*kki.*fphasei);
            bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*phasei.*fp;
            
            for j = 1:3
                for nii = 1: Nray
                    gradtempi = - ray(elem(:,i),nii);
                    gradtempi = gradtempi(:);
                    gradtempi = [real(gradtempi), imag(gradtempi)];
                    fphasei = gradtempi(:,1).*pxy(:,1) + gradtempi(:,2).*pxy(:,2);
                    ki = kk(elem(:,i),nii);
                    phasei = exp(1i*ki.*fphasei);
                    
                    for njj = 1: Nray
                        gradtempj = ray(elem(:,j),njj);
                        gradtempj = gradtempj(:);
                        gradtempj = [real(gradtempj), imag(gradtempj)];
                        fphasej = gradtempj(:,1).*pxy(:,1) + gradtempj(:,2).*pxy(:,2);
                        kj = kk(elem(:,j),njj);
                        phasej = exp(1i*kj.*fphasej);
                        exp_phase = phasei.*phasej;
                        
                        tempA1 = Dphi(:,1,i).*Dphi(:,1,j) + Dphi(:,2,i).*Dphi(:,2,j);
                        
                        tempA2 = 1i*ki*phi(p,i).*(gradtempi(:,1).*Dphi(:,1,j) + gradtempi(:,2).*Dphi(:,2,j))...
                            + 1i*kj*phi(p,j).*(gradtempj(:,1).*Dphi(:,1,i) + gradtempj(:,2).*Dphi(:,2,i));
                        
                        tempA3 = phi(p,i)*phi(p,j)*ki.*kj.*(gradtempi(:,1).*gradtempj(:,1)...
                            + gradtempi(:,2).*gradtempj(:,2));
                        
                        tempM = k2*phi(p,i)*phi(p,j);
                        
                        ii = (nii-1)*N + elem(:,i);
                        jj = (njj-1)*N + elem(:,j);
                        Aij = weight(p)*(tempA1 + tempA2 - tempA3 - tempM).*exp_phase.*area;
                        A1 = A1 + sparse(ii,jj,weight(p)*(tempA1).*exp_phase.*area,Ndof,Ndof);
                        A2 = A2 + sparse(ii,jj,weight(p)*(tempA2).*exp_phase.*area,Ndof,Ndof);
                        A3 = A3 + sparse(ii,jj,weight(p)*(tempA3).*exp_phase.*area,Ndof,Ndof);
                        M = M + sparse(ii,jj,weight(p)*(tempM).*exp_phase.*area,Ndof,Ndof);
                    end
                end
            end
        end
    end
    
    
    A = A1+A2-A3-M;
    
    B = A; B1 = A1; B2 = A2; B3 = A3; BM = M;
    
    
    %% %% 9th order Numerical Quadrature
    fquadorder = 9;     
    [lambda,weight] = quadpts(fquadorder);
    phi = lambda;           % the value of piesewise linear basis functions at numerical quadrature nodes
    nQuad = size(lambda,1);
    
    
    %% Compute geometric quantities and gradient of local basis
    [Dphi,area] = gradbasis(node,elem);
    reparea = repmat(area,Nray,1);
    
    
    %% Assemble the linear system Av = b
    
    
    
    A = sparse(Ndof,Ndof);
    A1 = A;
    A2 = A;
    A3 = A;
    M = A;
    
    bt = zeros(NT*Nray,3);
    b = zeros(Ndof,1);
    
    % fprintf('Assembling time:\n');
    % tic;
    for p = 1:nQuad
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        reppxy = repmat(pxy,Nray,1);
        fp = pde.f(reppxy);
        k2 = (omega./pde.speed(pxy)).^2;
        
        for i = 1:3
            gradtempi = - ray(elem(:,i),:);
            gradtempi = gradtempi(:);
            gradtempi = [real(gradtempi), imag(gradtempi)];
            fphasei = gradtempi(:,1).*reppxy(:,1) + gradtempi(:,2).*reppxy(:,2);
            kki = kk(elem(:,i),:);
            kki = kki(:);
            phasei = exp(1i*kki.*fphasei);
            bt(:,i) = bt(:,i) + weight(p)*phi(p,i)*phasei.*fp;
            
            for j = 1:3
                for nii = 1: Nray
                    gradtempi = - ray(elem(:,i),nii);
                    gradtempi = gradtempi(:);
                    gradtempi = [real(gradtempi), imag(gradtempi)];
                    fphasei = gradtempi(:,1).*pxy(:,1) + gradtempi(:,2).*pxy(:,2);
                    ki = kk(elem(:,i),nii);
                    phasei = exp(1i*ki.*fphasei);
                    
                    for njj = 1: Nray
                        gradtempj = ray(elem(:,j),njj);
                        gradtempj = gradtempj(:);
                        gradtempj = [real(gradtempj), imag(gradtempj)];
                        fphasej = gradtempj(:,1).*pxy(:,1) + gradtempj(:,2).*pxy(:,2);
                        kj = kk(elem(:,j),njj);
                        phasej = exp(1i*kj.*fphasej);
                        exp_phase = phasei.*phasej;
                        
                        tempA1 = Dphi(:,1,i).*Dphi(:,1,j) + Dphi(:,2,i).*Dphi(:,2,j);
                        
                        tempA2 = 1i*ki*phi(p,i).*(gradtempi(:,1).*Dphi(:,1,j) + gradtempi(:,2).*Dphi(:,2,j))...
                            + 1i*kj*phi(p,j).*(gradtempj(:,1).*Dphi(:,1,i) + gradtempj(:,2).*Dphi(:,2,i));
                        
                        tempA3 = phi(p,i)*phi(p,j)*ki.*kj.*(gradtempi(:,1).*gradtempj(:,1)...
                            + gradtempi(:,2).*gradtempj(:,2));
                        
                        tempM = k2*phi(p,i)*phi(p,j);
                        
                        ii = (nii-1)*N + elem(:,i);
                        jj = (njj-1)*N + elem(:,j);
                        Aij = weight(p)*(tempA1 + tempA2 - tempA3 - tempM).*exp_phase.*area;
                        A1 = A1 + sparse(ii,jj,weight(p)*(tempA1).*exp_phase.*area,Ndof,Ndof);
                        A2 = A2 + sparse(ii,jj,weight(p)*(tempA2).*exp_phase.*area,Ndof,Ndof);
                        A3 = A3 + sparse(ii,jj,weight(p)*(tempA3).*exp_phase.*area,Ndof,Ndof);
                        M = M + sparse(ii,jj,weight(p)*(tempM).*exp_phase.*area,Ndof,Ndof);
                    end
                end
            end
        end
    end
    
    
    A = A1+A2-A3-M;
    
    A = B-A; A1= B1-A1; A2 = B2-A2; A3 = B3-A3; M = BM-M;
    norm(A(:),inf)
    norm(A1(:),inf)
    norm(A2(:),inf)
    norm(A3(:),inf)
    norm(M(:),inf)
    
    
    
    
    
    
    
    
    
    
    
    toc;
end

% showrate(rec_omega(1:nn),rec_err1(1:nn)')
% showrate(rec_h(1:nn),rec_err1(1:nn))

fprintf(['\n' '-'*ones(1,80) '\n']);

fprintf('Error_1     ');  fprintf('  &  %.2e',rec_err1');  fprintf('\n\n');

