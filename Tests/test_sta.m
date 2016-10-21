numray3 = cell(N,1);
numray4 = cell(N,1);
dd = sqrt((node(:,1)-xs).^2 + (node(:,2)-ys).^2);
tic;
for i = 1:N
    temp = numray{i};
    numray3{i} = temp;
    numray4{i} = temp;
    
    if node(i,1) < 0 && node(i,2) > 0.2
        numray4{i} = [0, temp];
    end
    
    if dd(i) > 10*eps
        numray3{i} = [0, temp];
    end
end
toc;

[u3] = Ray_FEM_PML_22_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,numray3,fquadorder,plt);

[u4] = Ray_FEM_PML_22_PointSource(node,elem,omega,wpml,sigmaMax,xs,ys,speed,numray4,fquadorder,plt);




%% plot
a = 0.5;
u2 = reshape(u2,m,n);
u3 = reshape(u3,m,n);
u4 = reshape(u4,m,n);

dy = 0.8;
yn = round(dy/h) + 1;
xx = X(yn,:);
uu2 = u2(yn,:);
uu3 = u3(yn,:);
uu4 = u4(yn,:);

ryn = round(dy/fh) + 1;
rxx = -0.5:fh:0.5;
ruu = uf(ryn,:);

figure(21);
hold off;
plot(xx,real(uu2),'ro-');
hold on;
plot(rxx,real(ruu),'k');
xlabel('x');
ylabel('Wavefield');
legend('Ray-FEM solution','Reference solution','LOCATION','Best');
title('Wavefield at y = 0.3');

figure(22);
hold off;
plot(xx,real(uu3),'rs-');
hold on;
plot(rxx,real(ruu),'k');
xlabel('x');
ylabel('Wavefield');
legend('Ray-FEM solution with normal basis everywhere','Reference solution','LOCATION','Best');
title('Wavefield at y = 0.3');

figure(23);
hold off;
plot(xx,real(uu4),'r^-');
hold on;
plot(rxx,real(ruu),'k');
xlabel('x');
ylabel('Wavefield');
legend('Ray-FEM solution with normal basis at certain region','Reference solution','LOCATION','Best');
title('Wavefield at y = 0.3');



