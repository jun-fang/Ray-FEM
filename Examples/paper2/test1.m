%% set up
clear;
NPW = 4;                   % number of points per wavelength
test_num = 7;              % we test test_num examples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test 1: homogeneous with exact ray
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NPW = 4;                   
test1_hom_Exact_ray(NPW, test_num);

NPW = 6;                   
test1_hom_Exact_ray(NPW, test_num);

NPW = 8;                   
test1_hom_Exact_ray(NPW, test_num);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test 2: homogeneous with numerical ray
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NPW = 4;                  
test1_hom_Num_ray(NPW, test_num);

NPW = 6;                   
test1_hom_Num_ray(NPW, test_num);

NPW = 8;                  
test1_hom_Num_ray(NPW, test_num);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots 1: exact ray
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(11);
% 
% load('resutls_1_HomExRay_NPW_4.mat')
% show_convergence_rate(omegas(1:test_num), rel_l2_err(1:test_num),'omega','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');
% 
% load('resutls_1_HomExRay_NPW_6.mat')
% show_convergence_rate(omegas(1:test_num), rel_l2_err(1:test_num),'omega','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');
% 
% load('resutls_1_HomExRay_NPW_8.mat')
% show_convergence_rate(omegas(1:test_num), rel_l2_err(1:test_num),'omega','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots 2: numerical ray
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(12);
% 
% load('resutls_1_HomNumRay_NPW_4.mat')
% show_convergence_rate(omegas(1:test_num), rel_l2_err(1:test_num),'omega','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');
% 
% load('resutls_1_HomNumRay_NPW_6.mat')
% show_convergence_rate(omegas(1:test_num), rel_l2_err(1:test_num),'omega','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');
% 
% load('resutls_1_HomNumRay_NPW_8.mat')
% show_convergence_rate(omegas(1:test_num), rel_l2_err(1:test_num),'omega','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');


