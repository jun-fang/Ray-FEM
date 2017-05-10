%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test 1: homogeneous with exact ray
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
NPW = 4;                   % number of points per wavelength
test_num = 3;              % we test test_num examples
test1_hom_exact_ray(NPW, test_num);

clear;
NPW = 6;                   % number of points per wavelength
test_num = 3;              % we test test_num examples
test1_hom_exact_ray(NPW, test_num);

clear;
NPW = 8;                   % number of points per wavelength
test_num = 3;              % we test test_num examples
test1_hom_exact_ray(NPW, test_num);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test 2: homogeneous with numerical ray
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear;
% NPW = 4;                   % number of points per wavelength
% test_num = 6;              % we test test_num examples
% test2_hom_num_ray(NPW, test_num);

% clear;
% NPW = 6;                   % number of points per wavelength
% test_num = 6;              % we test test_num examples
% test2_hom_num_ray(NPW, test_num);
% 
% clear;
% NPW = 8;                   % number of points per wavelength
% test_num = 5;              % we test test_num examples
% test2_hom_num_ray(NPW, test_num);



