%% Test

NPW = 4;                   % number of points per wavelength
test_num = 7;              % we test test_num examples
test5_Lippmann_Schwinger(NPW, test_num);


NPW = 6;                   % number of points per wavelength
test_num = 7;              % we test test_num examples
test5_Lippmann_Schwinger(NPW, test_num);


NPW = 8;                   % number of points per wavelength
test_num = 7;              % we test test_num examples
test5_Lippmann_Schwinger(NPW, test_num);




%% Plot
% load('resutls_5_LipSch_NPW_4.mat');
% show_convergence_rate(high_omega(1:test_num),rel_l2_err(1:test_num),'omega','Rel L2 err');
% 
% load('resutls_5_LipSch_NPW_6.mat');
% show_convergence_rate(high_omega(1:test_num),rel_l2_err(1:test_num),'omega','Rel L2 err');
% 
% load('resutls_5_LipSch_NPW_8.mat');
% show_convergence_rate(high_omega(1:test_num),rel_l2_err(1:test_num),'omega','Rel L2 err');
