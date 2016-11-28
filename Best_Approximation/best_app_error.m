global omega a xs ys;
pde = Helmholtz_data1;
omega = 0*pi;
a = 1/2;
xs = 2; ys = 2;
fquadorder= 3;
plt = 0;
NPW = 6;
k = 1;

for oi = 1:1
    oi
    omega = omega + 10*pi;

h = 1/(NPW*round(omega/(2*pi)));
[node,elem] = squaremesh([-a,a,-a,a],h);
N = size(node,1); NT = size(elem,1); area = h*h/2;
ray = pde.ray(node);

rnode = node; relem = elem; fh = h;
for i = 1:k
    [rnode,relem] = uniformrefine(rnode,relem);
    fh = fh/2;
end

xnode = rnode(:,1); ynode = rnode(:,2);
rN = size(rnode,1);
Nray = size(ray,2);


% tic;
% ii = []; jj = []; ss = [];
% for nt = 1:NT
%     node_num = elem(nt,:);
%     node_nt = node(node_num,:);
%     rnt = elements_in_element(nt,node,elem,k);
%     rn = relem(rnt,:);
% %     idx = rn(:);
%     idx = unique(rn(:));
%     node1 = repmat(node_nt(1,:),size(idx));
%     node2 = repmat(node_nt(2,:),size(idx));
%     node3 = repmat(node_nt(3,:),size(idx));
%     lambda1 = tri_area(rnode(idx,:),node2,node3)/area;
%     lambda2 = tri_area(rnode(idx,:),node3,node1)/area;
%     lambda3 = tri_area(rnode(idx,:),node1,node2)/area;
%     mi = repmat(idx,1,3); mi = mi(:);
%     mj = repmat(node_num,size(idx)); mj = mj(:);
%     ms = [lambda1, lambda2, lambda3]; ms = ms(:);
%     idx = find(ms>0); mi = mi(idx); mj = mj(idx); ms = ms(idx);
%     ii = [ii;mi]; jj = [jj;mj]; ss = [ss;ms];  
% end
% 
% [~,I] = unique([ii,jj],'rows','first');
% ii = ii(I); jj = jj(I); ss = ss(I);
% toc;   
% 
% ii1 = 11;jj1=jj;ss1=ss;



tic;
ii = []; jj = []; ss = [];
for nt = 1:NT
    node_num = elem(nt,:);
    node_nt = node(node_num,:);
    xt = sort(node_nt(:,1)); yt = sort(node_nt(:,2));
    idx = find( (xnode<xt(3)+fh/2).*(xnode>xt(1)-fh/2)...
        .*(ynode<yt(3)+fh/2).*(ynode>yt(1)-fh/2) ); 
    node1 = repmat(node_nt(1,:),size(idx));
    node2 = repmat(node_nt(2,:),size(idx));
    node3 = repmat(node_nt(3,:),size(idx));
    lambda1 = tri_area(rnode(idx,:),node2,node3)/area;
    lambda2 = tri_area(rnode(idx,:),node3,node1)/area;
    lambda3 = tri_area(rnode(idx,:),node1,node2)/area;
    idxx = find((lambda1+lambda2+lambda3) <= (1+10*eps));
    idx = idx(idxx);
    mi = repmat(idx,1,3); mi = mi(:);
    mj = repmat(node_num,size(idx)); mj = mj(:);
    ms = [lambda1(idxx), lambda2(idxx), lambda3(idxx)]; ms = ms(:);
    idx = find(ms>0); mi = mi(idx); mj = mj(idx); ms = ms(idx);
    ii = [ii;mi]; jj = [jj;mj]; ss = [ss;ms];  
end

[~,I] = unique([ii,jj],'rows','first');
ii = ii(I); jj = jj(I); ss = ss(I);
toc;   

clear I;





% tic;
% P = sparse(rN,N);
% for nt = 1:NT
%     node_num = elem(nt,:);
%     node_nt = node(node_num,:);
%     node1 = repmat(node_nt(1,:),rN,1);
%     node2 = repmat(node_nt(2,:),rN,1);
%     node3 = repmat(node_nt(3,:),rN,1);
%     lambda1 = tri_area(rnode,node2,node3)/area;
%     lambda2 = tri_area(rnode,node3,node1)/area;
%     lambda3 = tri_area(rnode,node1,node2)/area;
%     idx = find(lambda1+lambda2+lambda3<=1+10*eps);
%     P(idx,node_num) = max(P(idx,node_num),[lambda1(idx), lambda2(idx), lambda3(idx)]);
% end
% toc;

tic;
NP = length(ii);
rows = []; cols = []; vals = [];
for nr = 1:Nray
    row = ii;
    col = jj + (nr-1)*NP;
    phase = real(ray(jj,nr)).*rnode(ii,1) + imag(ray(jj,nr)).*rnode(ii,2);
    val = ss.*exp(1i*omega*phase);
    rows = [rows; row];
    cols = [cols; col];
    vals = [vals; val];
end

A = sparse(rows,cols,vals,rN,N*Nray);
AT = transpose(A);

b = pde.ex_u(rnode);
u = (AT*A)\(AT*b);
toc;    

    %% Tirival way to build A
%     A = zeros(rN,N*Nray);
%     for nr = 1:Nray
%         for n = 1:N
%             xc = node(n,1); yc = node(n,2);
%             phi = nodal_basis(xc,yc,h,rnode);
%             phase = real(ray(n,nr))*rnode(:,1) + imag(ray(n,nr)).*rnode(:,2);
%             A(:,n+(nr-1)*N) = phi.*exp(1i*omega*phase);
%         end
%     end


[err1] = Ray_FEM_solution_error(node,elem,omega,pde,ray,u,fquadorder);

[~,~,v,~,~,err2] = Ray_FEM_IBC_1(node,elem,omega,pde,ray,fquadorder,plt);

err2/err1

end
