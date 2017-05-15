%% Test omega convergence rate for Lippmann-Schwinger equation

%% Set up
clear;
addpath(genpath('../../../ifem/'));
addpath('../../Functions/')
addpath('../../Methods/');
addpath('../../NMLA/');
addpath('../../Plots/');
% windows path
addpath('C:\Users\Jun Fang\Documents\MATLAB\Solutions_Lippmann_Schwinger');
% linux path
% addpath('/home/jun/Documents/MATLAB/Solutions_Lippmann_Schwinger');

test_num = 7;              % we test test_num examples


NPW = 4;                   % number of points per wavelength
test2_Lippmann_Schwinger(NPW, test_num);


NPW = 6;                   % number of points per wavelength
test2_Lippmann_Schwinger(NPW, test_num);


NPW = 8;                   % number of points per wavelength
test2_Lippmann_Schwinger(NPW, test_num);




%% Plot
% load('resutls_2_LipSch_NPW_4.mat');
% show_convergence_rate(high_omega(1:test_num),rel_l2_err(1:test_num),'omega','Rel L2 err');
% 
% load('resutls_2_LipSch_NPW_6.mat');
% show_convergence_rate(high_omega(1:test_num),rel_l2_err(1:test_num),'omega','Rel L2 err');
% 
% load('resutls_2_LipSch_NPW_8.mat');
% show_convergence_rate(high_omega(1:test_num),rel_l2_err(1:test_num),'omega','Rel L2 err');

