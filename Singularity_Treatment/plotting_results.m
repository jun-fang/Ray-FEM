
load('resutls_1_HomExRay_NPW_4.mat')

test_num = 9

figure(22);
show_convergence_rate(omegas(1:test_num), rel_l2_err,'omega','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');

 
load('resutls_1_HomExRay_NPW_6.mat')

test_num = 8
 
 
show_convergence_rate(omegas(1:test_num), rel_l2_err,'omega','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');

  
load('resutls_1_HomExRay_NPW_8.mat')

test_num  =7

show_convergence_rate(omegas(1:test_num), rel_l2_err,'omega','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');



load('resutls_2_HomNumRay_NPW_4.mat')

figure(1); clf(); 
show_convergence_rate(high_omega(1:test_num),rel_l2_err(1:test_num),'omega','Rel L2 ');
 hold on; 
 

load('resutls_2_HomNumRay_NPW_6.mat')
show_convergence_rate(high_omega(1:test_num),rel_l2_err(1:test_num),'omega','Rel L2 ');


load('resutls_2_HomNumRay_NPW_8.mat')
show_convergence_rate(high_omega(1:test_num),rel_l2_err(1:test_num),'omega','Rel L2 ');


