function [ub, ub_g1, ub_g2] = Babich_expansion(xs,ys,pxy,omega,option,Bx,By,phase,amplitude)
%% Construct Babich expansion and its gradient


% homogeneous case
if ~iscell(option) && strcmp(option, 'homogeneous')
    x = (pxy(:,1)-xs);  y = (pxy(:,2)-ys);
    r = sqrt(x.^2 + y.^2);
    
    ub = 1i/4*besselh(0,1,omega*r);                  % Babich expression
    ub_g1 = -1i/4*omega*besselh(1,1,omega*r)./r.*x;  % partial derivative wrt x
    ub_g2 = -1i/4*omega*besselh(1,1,omega*r)./r.*y;  % partial derivative wrt y
end

% gravity case
if iscell(option) && strcmp(option{1}, 'gravity')
    alpha = option{2};   E = option{3};
    trg = pxy';  src = [xs;ys];
    ub = lhelmfs(trg,src,alpha,E);
    
    trg = repmat([xs;ys], 1, size(pxy,1)); src = pxy';
    [~, ub_g1, ub_g2] = lhelmfs(trg,src,alpha,E,1);
    ub = ub(:);  ub_g1 = ub_g1(:);  ub_g2 = ub_g2(:);
end

% Babich case
if  iscell(option) && strcmp(option{1}, 'Babich')
    
    Bu = amplitude;  %[DD1(:), DD1x(:), DD1y(:), DD2(:), DD2x(:), DD2y(:)];
    Buint = interpolation2(Bx, By, Bu, pxy);
    amp1 = Buint(:,1);  amp1x = Buint(:,2);  amp1y = Buint(:,3);
    amp2 = Buint(:,4);  amp2x = Buint(:,5);  amp2y = Buint(:,6);
    
    if strcmp(option{3}, 'numerical_phase')
        Bu = phase; %[ttao(:), ttaox(:), ttaoy(:)];
        Buint = interpolation2(Bx, By, Bu, pxy);
        pha = Buint(:,1);  phax = Buint(:,2);  phay = Buint(:,3);
    end
    
    if strcmp(option{3}, 'exact_phase')
        if strcmp(option{2}, 'CGV')
            [~, pha, phax, phay] = eikonal_cgv(1, [0.1, -0.2], [0,0], pxy);
        end
        
        if strcmp(option{2}, 'Homo')
            x = (pxy(:,1)-xs);  y = (pxy(:,2)-ys);
            pha = sqrt(x.^2 + y.^2);
            phax = x./pha;  phay = y./pha;
        end
    end
    
    
    c1 = 1i*(sqrt(pi)/2);
    f1 = c1*besselh(0,1,omega*pha);
    f1x = - c1*besselh(1,1,omega*pha)*omega.*phax;
    f1y = - c1*besselh(1,1,omega*pha)*omega.*phay;
    
    G1 = f1.*amp1;
    G1x = f1x.*amp1 + f1.*amp1x;
    G1y = f1y.*amp1 + f1.*amp1y;
    
    c2 = 1i*(sqrt(pi)/2)*exp(1i*pi);
    f2 = c2*(2*pha/omega).*besselh(1,1,omega*pha);
    temp = c2* ( 4/omega*besselh(1,1,omega*pha) - 2*pha.*besselh(2,1,omega*pha) );
    f2x = temp.*phax;
    f2y = temp.*phay;
    
    G2 = f2.*amp2;
    G2x = f2x.*amp2 + f2.*amp2x;
    G2y = f2y.*amp2 + f2.*amp2y;
    
    ub = G1 + G2;
    ub_g1 = G1x + G2x;
    ub_g2 = G1y + G2y;
    
end

