%% Test omega convergence rate for Lippmann-Schwinger equation

%% Set up
clear;
addpath(genpath('../../../ifem/'));
addpath('../../Functions/')
addpath('../../Methods/');
addpath('../../NMLA/');
addpath('../../Plots/');
addpath('../../Helmholtz_data/');

test_num = 7;              % we test test_num examples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test 1: constant gradient of velocity with exact ray
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NPW = 4;                   % number of points per wavelength
% test3_CGV_Exact_ray(NPW, test_num);
% 
% 
% NPW = 6;                   
% test3_CGV_Exact_ray(NPW, test_num);
% 
% 
% NPW = 8;                   
% test3_CGV_Exact_ray(NPW, test_num);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test 2: constant gradient of velocity with numerical ray
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NPW = 4;                   
% test3_CGV_Num_ray(NPW, test_num);
% 
% 
% NPW = 6;                   
% test3_CGV_Num_ray(NPW, test_num);
% 
% 
% NPW = 8;                   
% test3_CGV_Num_ray(NPW, test_num);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots 1: exact ray
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(31);
% 
% load('results_3_CGV_ExRay_NPW_4.mat')
% show_convergence_rate(high_omega(1:test_num), rel_l2_err(1:test_num),'omega','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');
% 
% load('results_3_CGV_ExRay_NPW_6.mat')
% show_convergence_rate(high_omega(1:test_num), rel_l2_err(1:test_num),'omega','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');
% 
% load('results_3_CGV_ExRay_NPW_8.mat')
% show_convergence_rate(high_omega(1:test_num), rel_l2_err(1:test_num),'omega','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plots 2: numerical ray
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure(32);
% 
% load('results_3_CGV_NumRay_NPW_4.mat')
% show_convergence_rate(high_omega(1:test_num), rel_l2_err(1:test_num),'omega','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');
% 
% load('results_3_CGV_NumRay_NPW_6.mat')
% show_convergence_rate(high_omega(1:test_num), rel_l2_err(1:test_num),'omega','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');
% 
% load('results_3_CGV_NumRay_NPW_8.mat')
% show_convergence_rate(high_omega(1:test_num), rel_l2_err(1:test_num),'omega','||u - u_h||_{L^2(\Omega)}/||u||_{L^2(\Omega)}');

