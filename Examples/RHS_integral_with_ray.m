function b = RHS_integral_with_ray(xc,yc,h,rh,omega,ray,source,fquadorder)

[node,elem] = squaremesh([xc-h,xc+h,yc-h,yc+h],rh);


%% Numerical Quadrature
[lambda,weight] = quadpts(fquadorder);
nQuad = size(lambda,1);


%% Compute geometric quantities and gradient of local basis
[~,area] = gradbasis(node,elem);

%% Assemble right-hand side

bt = zeros(size(elem,1),1);       % the right hand side

for p = 1:nQuad
    % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:);
    
    % we suppose that the source is well inside the physical domain
    fp = source(pxy);
    fphase = real(ray)*pxy(:,1) + imag(ray)*pxy(:,2);
    exp_phase = exp(-1i*omega*fphase);
%     fp = source(pxy).*( pxy(:,1) < xmax - wpml ).*( pxy(:,1) > xmin + wpml )...
%         .*( pxy(:,2) < ymax - wpml ).*( pxy(:,2) > ymin + wpml ); 
    bt = bt + weight(p)*nodal_basis(xc,yc,h,pxy).*exp_phase.*fp.*area;
end

b = sum(bt);

