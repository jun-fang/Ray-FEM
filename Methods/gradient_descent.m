function [opt_angs,opt_Bs,opt_res,niter,org_res] = gradient_descent(node,u,omega,x0,y0,c0,angs,Bs)

% termination tolerance
tol = 1e-6;

% maximum number of allowed iterations
maxiter = 20;          
   

% step size ( 0.33 causes instability, 0.2 quite accurate)
alpha = 1;
beta = 0.5;
c1 = 1e-4;
c2 = 0.9;

% initialize gradient norm, optimization vector, iteration counter, perturbation
res = residual(node,u,omega,x0,y0,c0,angs,Bs);
org_res = res;
niter = 0; 
d_res = 1; 

% residual(node,u,omega,speed,x0,y0,pi/4,1)

% gradient descent algorithm:
while and(and(res>=tol,d_res>=tol), niter <= maxiter)
    
    niter = niter + 1;
    
    % calculate gradient:
    [gd1, gd2] = grad(node,u,omega,x0,y0,c0,angs,Bs);
    new_angs = angs - alpha*gd1;
    new_Bs = Bs - alpha*gd2;  
    new_res = residual(node,u,omega,x0,y0,c0,new_angs,new_Bs);
%     [new_gd1, new_gd2] = grad(node,u,omega,x0,y0,c0,new_angs,new_Bs);

    %% Wolfe condition: Armijo/backtracking line search + curvature condition 
    while new_res > res - c1*alpha*(norm(gd1,2)^2 + norm(gd2,2)^2)
%             abs(gd1.*new_gd1 + gd2.*new_gd2) < c2*(norm(gd1,2)^2 + norm(gd2,2)^2) )
        alpha = beta*alpha;
        new_angs = angs - alpha*gd1;
        new_Bs = Bs - alpha*gd2;  
        new_res = residual(node,u,omega,x0,y0,c0,new_angs,new_Bs);
%         [new_gd1, new_gd2] = grad(node,u,omega,x0,y0,c0,new_angs,new_Bs);
    end
    
%     alpha = alpha
%     new_res = new_res
    d_res = (res - new_res)/res;
    angs = new_angs;
    Bs = new_Bs;
    res = new_res;
    
    
    
end

niter = niter - 1;
opt_angs = angs;
opt_Bs = Bs;
opt_res = res;

     


%% u difference
    function du = u_diff(node,u,omega,x0,y0,c0,x,y,angs,Bs)
        N = length(angs);
        temp = zeros(size(x));
        for n = 1:N
            dot_temp = cos(angs(n))*(x-x0) + sin(angs(n))*(y-y0);
            temp = temp + Bs(n)*exp(1i*omega/c0*dot_temp);
        end
        idx = xy_to_index(node,x,y);
        du = u(idx) - temp;
    end

%% residual
    function res = residual(node,u,omega,x0,y0,c0,angs,Bs)
        h = node(2,2) - node(1,2);
        
%         rd = 9/omega^(3/4);
%         M = 32;
%         angl = linspace(0,2*pi,M+1) ; 
%         ang = angl(1:M) ;
%         x = x0 + rd*cos(ang');
%         y = y0 + rd*sin(ang');
%         [~, x, y] = xy_to_index(node,x,y);
        r1 = 6/omega^(3/4);
        r2 = 8/omega^(3/4);
        theta1 = 7/6*pi;
        theta2 = 4/3*pi;
        [rr,theta] = meshgrid(r1:2*h:r2,theta1:pi/16:theta2);
        x = x0 + rr.*cos(theta);
        y = y0 + rr.*sin(theta);
        x = x(:);
        y = y(:);
        
        ud = u_diff(node,u,omega,x0,y0,c0,x,y,angs,Bs);
        res = norm(ud,2);
    end

%% gradient
    function [gd1, gd2] = grad(node,u,omega,x0,y0,c0,angs,Bs)
        % gd1: partial derivatie wrt angle theta_n
        % gd2: partial derivatie wrt amplitude B_n
        N = length(angs);
        gd1 = zeros(1,N);   gd2 = zeros(1,N);
        
        h = node(2,2) - node(1,2);
        ikc0 = 1i*omega/c0;
        
%         rd = 9/omega^(3/4);
%         M = 32;
%         angl = linspace(0,2*pi,M+1) ; 
%         ang = angl(1:M) ;
%         x = x0 + rd*cos(ang');
%         y = y0 + rd*sin(ang');
        
        r1 = 6/omega^(3/4);
        r2 = 8/omega^(3/4);
        theta1 = 7/6*pi;
        theta2 = 4/3*pi;
        [rr,theta] = meshgrid(r1:2*h:r2,theta1:pi/16:theta2);
        x = x0 + rr.*cos(theta);
        y = y0 + rr.*sin(theta);
        x = x(:);
        y = y(:);
        
        [~, x, y] = xy_to_index(node,x,y);
        
        ud = u_diff(node,u,omega,x0,y0,c0,x,y,angs,Bs);
        
        for n = 1:N
            temp1 =  -sin(angs(n))*(x-x0) + cos(angs(n))*(y-y0) ;
            temp2 = exp( ikc0*( cos(angs(n))*(x-x0) + sin(angs(n))*(y-y0) ) );
            
            gd1(n) = 2*sum(real( -Bs(n)*temp2*ikc0.*temp1.*conj(ud) ));
            gd2(n) = -2*sum( conj(temp2.*ud) );

        end
    end

end


