test_MPM10_2_point_sources_two_directions;
test_MPM10_3_point_sources_two_directions;
test_MPM10_4point_sources_two_directions;




% rec_omega = [25 50 100 200 400 800 1600]*pi;
% rec_err1 = zeros(length(rec_omega),1);
% rec_err2 = rec_err1;
% rec_err3 = rec_err1;
% rec_err0 = rec_err1;
% global omega;
% global xs;
% global ys;
% global a;
% 
% pde = Helmholtz_data1_0;
% plt = 0;                           % plot solution or not
% fquadorder = 3;                    % numerical quadrature order
% solver = 'DIR';
% 
% xs = 10;
% ys = 10;
% a = 1/2;
% 
% nn = 2;
% 
% for j = 1:nn%length(rec_omega)
%    
%     low_omega = sqrt(rec_omega(j));
%     omega = rec_omega(j);
%     NPW = 8;
%     h = 1/round((NPW*omega)/(2*pi));   % mesh size
%     
%     omega = low_omega;
%     
%     wl = 2*pi/omega;    % wavelength
%     rh = wl/NPW;
%     
%     a = 1/2;
%     a = a + h*(round(wl/h) + 1);
%     a
%     1/h
%     [lnode,lelem] = squaremesh([-a,a,-a,a],h);
%     [lu,~,~,rel_L2_err] = Standard_FEM_IBC(lnode,lelem,omega,pde,fquadorder,solver,plt);
%     eu = pde.ex_u(lnode);
% 
% %     rec_err0(j) = rel_L2_err; %norm(lu-eu,inf);
%     rec_err0(j) = norm(lu-eu,inf);
% 
% end
% 
% fprintf(['\n' '-'*ones(1,80) '\n']);
% fprintf('Error_0     ');
% fprintf('  &  %.2e',rec_err0');
% fprintf('\n\n');
% 
% 
% showrate(rec_omega(1:nn), rec_err0(1:nn)')
