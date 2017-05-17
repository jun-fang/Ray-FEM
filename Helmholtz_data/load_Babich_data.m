function [Bh0,Bx0,By0,D1,D2,tao,tao2x,tao2y] = load_Babich_data(omega, option)

if strcmp(option, 'CGV')
    switch round(omega/(pi))
        case 100
            load('Babich_CGV_40.mat');
        case 150
            load('Babich_CGV_60.mat');
        case 200
            load('Babich_CGV_80.mat');
        case 300
            load('Babich_CGV_120.mat');
        case 500
            load('Babich_CGV_200.mat');
        case 750
            load('Babich_CGV_300.mat');
        case 1000
            load('Babich_CGV_300.mat');
    end
end

if strcmp(option, 'Homo')
    switch round(omega/(pi))
        case 250
            load('Babich_Homo_25.mat');
        case 400
            load('Babich_Homo_40.mat');
        case 600
            load('Babich_Homo_60.mat');
        case 1000
            load('Babich_Homo_100.mat');
        case 1500
            load('Babich_Homo_150.mat');
    end
end