function [A,b] = modify_linear_system_with_boundary_conditons(A,b,xc,yc,node,elem,omega,speed,ray,fquadorder,opt)

N = size(node,1);       % number of grid points
Nray = size(ray,2);     % number of rays crossing at each grid node
Ndof = N*Nray;          % degree of freedom

c = speed(node);    % medium speed
k = omega./c;           % wavenumber
kk = repmat(k,1,Nray);

xmax = max(node(:,1));
xmin = min(node(:,1));
ymax = max(node(:,2));
ymin = min(node(:,2));

%% Boundaries
[~,bdEdge,~] = findboundary(elem);
bdxy = node(bdEdge(:,1),:);
bx = bdxy(:,1); by = bdxy(:,2);
idx = ( bx>xmin+10*eps ).*( bx<xmax-10*eps )...
    .*( by>ymin+10*eps ).*( by<ymax-10*eps );
bdEdge = bdEdge(idx>0,:);
bdxy = node(bdEdge(:,1),:);
x1 = max(bdxy(:,1));
x2 = min(bdxy(:,1));
y1 = max(bdxy(:,2));
y2 = min(bdxy(:,2));

Ne = size(bdEdge,1);      % number of boundary edges
el = sqrt(sum((node(bdEdge(:,1),:) - node(bdEdge(:,2),:)).^2,2));
repel = repmat(el,Nray,1);

[lambdagN,weightgN] = quadpts1(fquadorder);
phigN = lambdagN;
nQuadgN = size(lambdagN,1);


%% Neumann Boundary condition
if strcmp(opt,'Neu')
    ge = zeros(Ne*Nray,2);
    for pp = 1:nQuadgN
        ppxy = lambdagN(pp,1)*node(bdEdge(:,1),:) ...
            + lambdagN(pp,2)*node(bdEdge(:,2),:);
        kxy = omega./speed(ppxy);
        reppxy = repmat(ppxy,Nray,1);
        gp = g_Neu(omega,reppxy,xc,yc,x1,x2,y1,y2);
        for i = 1:2
            gradtempi = - ray(bdEdge(:,i),:);
            gradtempi = gradtempi(:);
            gradtempi = [real(gradtempi), imag(gradtempi)];
            fphasei = gradtempi(:,1).*reppxy(:,1) + gradtempi(:,2).*reppxy(:,2);
            kki = kk(bdEdge(:,i),:);
            kki = kki(:);
            phasei = exp(1i*kki.*fphasei);
            ge(:,i) = ge(:,i) + weightgN(pp)*phigN(pp,i)*phasei.*gp;
        end
    end
    ge = ge.*repmat(repel,1,2);
    for ni = 1:Nray
        ii = (ni-1)*Ne+1:1:ni*Ne;
        jj = (ni-1)*N+1:1:ni*N;
        geii = ge(ii,:);
        b(jj) = b(jj) + accumarray(bdEdge(:), geii(:),[N 1]);
    end
end








%% Modify the linear system Av = b with Impedance Boundary Condition
if strcmp(opt, 'IBC')
    rows = zeros(nQuadgN*2*2*Nray*Nray*Ne, 1);
    cols = rows;
    vals = rows;
    inds = 1:Ne;
    
    ge = zeros(Ne*Nray,2);
    
    for pp = 1:nQuadgN
        ppxy = lambdagN(pp,1)*node(bdEdge(:,1),:) ...
            + lambdagN(pp,2)*node(bdEdge(:,2),:);
        kxy = omega./speed(ppxy);
        reppxy = repmat(ppxy,Nray,1);
        gp = g_IBC(omega,reppxy,xc,yc,x1,x2,y1,y2);
        
        for i = 1:2
            gradtempi = - ray(bdEdge(:,i),:);
            gradtempi = gradtempi(:);
            gradtempi = [real(gradtempi), imag(gradtempi)];
            fphasei = gradtempi(:,1).*reppxy(:,1) + gradtempi(:,2).*reppxy(:,2);
            kki = kk(bdEdge(:,i),:);
            kki = kki(:);
            phasei = exp(1i*kki.*fphasei);
            ge(:,i) = ge(:,i) + weightgN(pp)*phigN(pp,i)*phasei.*gp;
            
            for j = 1:2
                for nii = 1:Nray
                    gradtempi = - ray(bdEdge(:,i),nii);
                    gradtempi = gradtempi(:);
                    gradtempi = [real(gradtempi), imag(gradtempi)];
                    fphasei = gradtempi(:,1).*ppxy(:,1) + gradtempi(:,2).*ppxy(:,2);
                    ki = kk(bdEdge(:,i),nii);
                    phasei = exp(1i*ki.*fphasei);
                    
                    for njj = 1: Nray
                        gradtempj = ray(bdEdge(:,j),njj);
                        gradtempj = gradtempj(:);
                        gradtempj = [real(gradtempj), imag(gradtempj)];
                        fphasej = gradtempj(:,1).*ppxy(:,1) + gradtempj(:,2).*ppxy(:,2);
                        kj = kk(bdEdge(:,j),njj);
                        phasej = exp(1i*kj.*fphasej);
                        exp_phase = phasei.*phasej;
                        
                        rows(inds) = (nii-1)*N + bdEdge(:,i);
                        cols(inds) = (njj-1)*N + bdEdge(:,j);
                        vals(inds) = 1i*kxy*weightgN(pp)*phigN(pp,i)*phigN(pp,j).*exp_phase.*el;
                        inds = inds + Ne;
                    end
                end
            end
        end
    end
    
    A = A + sparse(rows,cols,vals,Ndof,Ndof);
    
    ge = ge.*repmat(repel,1,2);
    for ni = 1:Nray
        ii = (ni-1)*Ne+1:1:ni*Ne;
        jj = (ni-1)*N+1:1:ni*N;
        geii = ge(ii,:);
        b(jj) = b(jj) + accumarray(bdEdge(:), geii(:),[N 1]);
    end
    
end


