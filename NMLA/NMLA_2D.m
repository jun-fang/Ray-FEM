function [angs] = NMLA_2D(x0,y0,c0,omega,Rest,node,elem,u,ux,uy,pde,pct,Nray,opt,plt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NMLA: finding ray directions at observation point (x0,y0) in 2D case
%  See details in 'Numerical MicroLocal Analysis Revisited' by Jean-David
% Benamou, Francis Collino, SimonMarmorat : https://hal.inria.fr/inria-00558881/document
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  INPUT:
% 
%     (x0,y0):     Observation point; 
%     c0:          Medium speed at (x0,y0)
%     omega:       Angular frequency
% 
%     Rest:        The estimate of the distance between the source point and
%                the observation point. 
%
%     node:        N x 2 matrix that contains the physical position of each node
%                node(:,1) provides the x coordinates, node(:,2) provides the y coordinates      
%
%     elem:        NT x 3 matrix that contains the indices of the nodes for each
%                edged of each triangle 
%
%     u,ux,uy:     N x 1 vectors. Solution data u, \partial_x u, \partial_y u
%                at each grid point, these could be numerically computed data
%
%     pde:         Structure of functions containing the pde information 
%                like the exact solution u(x), gradient Du(x), etc.
% 
%     pct:         The percentage of the highest peak to claim a peak.
%                For example, if pct = 1/2, then all the peaks should be
%                above 1/2 of the highest peak. Otherwise, we don't count
%                it as a peak.
% 
%     Nray:        Number of ray directions (may not be needed)
%                Nray = 0 means we don't know how many ray directions.
%                If Nray is given, we can capture ray directions more effeciently. 
% 
%     opt:         If opt = 'num', we use numerically computed data u,ux,uy 
%                to interpolate these data on the observation circle; 
%                  If opt = 'ex', we use the exact solution u(x) and exact
%                gradient Du(x) on the observation circle.
%
%     plt:         If plt = 1, then show the plot; if plt = 0, don show it. 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  OUTPUT:
%
%     angs:        Dominant ray direction angles
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameters
gamma = 0.5 ;               %% parameter in impedance quantity
largd = 3.5 ;               %% Parameters in the Gaussian function

%% Truncation level and number of samples
p = [1,1,-2.5-0.93*(omega*Rest)^0.763];
rt = roots(p);
pidx = find(rt>0);
pr = rt(pidx(1))/omega;     %% real radius of observation circle
kr0 = pr*omega/c0;          %% k*r0
rlel = round(kr0 + (kr0)^(1/3) -2.5);    %% truncation level to obtain needed precision
rlel = max(1,rlel) ;
nb_theta = 2*(5*rlel)+1 ;                %% number of samples on the observation circle

%% Angle discretizaion on the circle
angl = linspace(0,2*pi,nb_theta+1) ;  
ang=angl(1:nb_theta) ;
X = x0+ pr*cos(ang) ;
Y = y0 + pr*sin(ang) ;

%% Data on the circle
if (strcmp(opt,'num'))     %% use numerically computed data
    
    xmin = node(1,1);    xmax = node(end,1);
    ymin = node(1,2);    ymax = node(end,2);
    dx = min(x0-xmin, xmax-x0);
    dy = min(y0-ymin, ymax-y0);
    if dx < pr || dy < pr
        fprintf('Observation circle out of the computational domain: pr = %d\n\n', pr);
        angl = linspace(0,2*pi,8);
        angs = angl(1:7);
        return;
    end
      
    Field = interpolation(node,elem,[X',Y'],u);
    DUx = interpolation(node,elem,[X',Y'],ux);
    DUy = interpolation(node,elem,[X',Y'],uy);
end

if (strcmp(opt,'ex'))      %% use the exact data
    xynode = [X',Y'];
    Field = pde.ex_u(xynode);
    DU = pde.Du(xynode);
    DUx = DU(:,1);
    DUy = DU(:,2);
end

DField = DUx.*cos(ang') + DUy.*sin(ang');
U = gamma*DField/(1i*omega/c0) + Field;
U = transpose(U);

%% Filtering 
fu = fft(U) ;
[fbeta] = BGFiltrage(fu,kr0,gamma,rlel,largd,nb_theta) ;
beta = ifft(fbeta) ;       %% amplitudes in the direction given by ang

%% Plot
if plt
    figure(1);
    plot(ang,abs(beta),'r-');
end

%% Find dominant peaks and ray direction angles
[mm,ii] = max(abs(beta));  %% find max and significant angle
[pks,locs] = findpeaks(abs(beta),'MinPeakHeight',pct*mm);
if Nray == 1               %% given one ray direciton
    locs = ii;
    angs = ang(locs);
    Bs = mm;
elseif Nray                %% given multiple ray directions
    [~, high_idx] = sort(-pks);
    idx = high_idx(1:Nray);
    locs = locs(idx);
    Bs = pks(idx);
    [locs,idx] = sort(locs);
    Bs = Bs(idx);
    angs = ang(locs);
else                       %% number of ray directions is not given
    angs = ang(locs);
    Bs = pks;
end

%% Post processing 
% lb = [];
% ub = [];
% options = optimset('Display', 'off');
% x = lsqnonlin(@(x) residual_fun(x), [angs;Bs], lb, ub, options);
% NN = round(length(x)/2);
% x = x(1:NN);
% angs = sort(x);
% 
% 
% function res = residual_fun(x)
% NN = round(length(x)/2);
% res = zeros(nb_theta,1);
% for n = 1:NN
%     res = res + x(NN + n)*exp(1i*kr0*( cos(ang')*cos(x(n)) + sin(ang')*sin(x(n))));
% end
% res = Field - res;
% end

end
