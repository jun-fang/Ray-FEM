function [angs, r, Bs] = NMLA(x0,y0,c0,omega,Rest,node,elem,u,ux,uy,pde,pct,Nray,data,opt,plt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NMLA 2nd order correction: finding ray directions at observation point (x0,y0) in 2D case
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
%     data:        If data = 'num', we use numerically computed data u,ux,uy
%                to interpolate these data on the observation circle;
%                  If data = 'ex', we use the exact solution u(x) and exact
%                gradient Du(x) on the observation circle.
%
%     opt:         If opt = 1, we do sencord order correction;
%                  If opt = 0, we do the normal one.
%
%     plt:         If plt = 1, then show the plot; if plt = 0, don show it.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  OUTPUT:
%
%     angs:        Dominant ray direction angles
%     r:           Radius of the observation circle
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parameters
imp = 0.5 ;                %% parameter in impedance quantity
gau = 3.5 ;                %% Parameter in Gaussian function


%% Truncation level and number of samples
r = NMLA_radius(omega/c0,Rest);
kr = r*omega/c0;                      %% k*r
L = round(kr + (kr)^(1/3) -2.5);      %% truncation level to obtain needed precision
L = max(1,L) ;
M = 2*(8*L)+1;                        %% number of samples on the observation circle


%% Angle discretizaion on the circle
angl = linspace(0,2*pi,M+1) ;
ang=angl(1:M) ;
X = x0+ r*cos(ang) ;
Y = y0 + r*sin(ang) ;


%% Data on the circle
if (strcmp(data,'num'))      %% use numerically computed data
    xmin = node(1,1);    xmax = node(end,1);
    ymin = node(1,2);    ymax = node(end,2);
    
    dx = min(x0-xmin, xmax-x0);
    dy = min(y0-ymin, ymax-y0);
    if dx < r || dy < r
        fprintf('Observation circle out of the computational domain: r = %d\n\n', r);
        angl = linspace(0,2*pi,8);
        angs = angl(1:7);
        return;
    end
    
    Field = interpolation(node,elem,[X',Y'],u);
    DUx = interpolation(node,elem,[X',Y'],ux);
    DUy = interpolation(node,elem,[X',Y'],uy);
    
end

if (strcmp(data,'ex'))        %% use the exact data
    xynode = [X',Y'];
    Field = pde.u_ex(xynode);
    DU = pde.Du(xynode);
    DUx = DU(:,1);
    DUy = DU(:,2);
end

DField = DUx.*cos(ang') + DUy.*sin(ang');
U = imp*DField/(1i*omega/c0) + Field;
U = transpose(U);


%% Filtering
fu = fft(U) ;
[fbeta] = BGFiltrage(fu,kr,imp,L,gau,M) ;
beta = ifft(fbeta) ;


%% Plot
if plt
    figure(2);
    plot(ang,abs(beta),'r-');
end

%% Extend beta to include the endpoints
ex_beta = [beta,beta(1:2)];

%% Find dominant peaks and ray direction angles
[mm,ii] = max(abs(ex_beta));    %% find max and significant angle
% [~,locs] = findpeaks(abs(ex_beta),'MinPeakDistance',M/32,'MinPeakHeight',pct*mm);
[~,locs] = findpeaks(abs(ex_beta),'MinPeakDistance',M/24,'MinPeakHeight',pct*mm);
locs = (locs>M).*(locs-M) + (locs<M+1).*locs;

pks = beta(locs);

if Nray == 1                 %% given one ray direction
    locs = ii;
    angs = ang(locs);
    Bs = mm;
elseif Nray > 1                 %% given multiple ray directions
    [~, high_idx] = sort(-pks);
    Nray = min(Nray, length(high_idx));
    idx = high_idx(1:Nray);
    locs = locs(idx);
    Bs = pks(idx);
    [locs,idx] = sort(locs);
    Bs = Bs(idx);
    angs = ang(locs);
elseif Nray == 0                        %% number of ray directions is not given
    angs = ang(locs);
    Bs = pks;
    if length(angs)>3   %% one can manually set how many rays at most
        [~, high_idx] = sort(-pks);
        idx = high_idx(1:3);
        locs = locs(idx);
        Bs = pks(idx);
        [locs,idx] = sort(locs);
        Bs = Bs(idx);
        angs = ang(locs);
    end
end


%% Second order correction: for one ray case
if (opt)
    ang_est = angs(1);
    %     rts = roots([1, 4*omega*pi*Rest/M, -2*omega*pi*Rest]);
    %     rt1 = max(rts);
    L = min(L,round(sqrt(2*pi*omega*Rest)/4));        % such that the imaginary part would not excced \pi
    lbeta = zeros(2*L+1,1);
    ii = 1:L;
    jj = -L:1:-1;
    lbeta(L+1) = fbeta(1);
    lbeta(L+1+ii) = fbeta(ii+1).*exp(1i*ii*ang_est);
    lbeta(1:L) = fbeta(M+1+jj).*exp(1i*jj*ang_est);
    yy = imag(log(lbeta/lbeta(L+1)));
    xx = -L:1:L;
    %     plot(xx,yy);
    pf = polyfit(xx',yy,2);          %% complexity O(L)
    angs(1) = angs(1) - pf(2);
end



end
