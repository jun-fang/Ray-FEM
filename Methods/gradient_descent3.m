function [opt_angs,opt_Bs,opt_res,niter,org_res] = gradient_descent2(node,u,omega,x0,y0,c0,angs,Bs)

% termination tolerance
tol = 1e-6;

% maximum number of allowed iterations
maxiter = 20;          
   

% step size ( 0.33 causes instability, 0.2 quite accurate)
alpha = 0.2;
beta = 0.5;

% initialize gradient norm, optimization vector, iteration counter, perturbation
res = residual(node,u,omega,x0,y0,c0,angs,Bs)
org_res = res;
niter = 0; 
d_res = 1; 

% residual(node,u,omega,speed,x0,y0,pi/4,1)

% gradient descent algorithm:
while and(and(res>=tol,d_res>=tol), niter <= maxiter)
    
    niter = niter + 1
    
    % calculate gradient:
    [gd1, gd2] = grad(node,u,omega,x0,y0,c0,angs,Bs)
%     gd1 = 10*pi*gd1;
    t = 1;
    new_angs = angs - t*gd1;
    new_Bs = Bs - t*gd2;  
    new_res = residual(node,u,omega,x0,y0,c0,new_angs,new_Bs);

    %% Armijo/backtracking line search 
    while new_res > res - alpha*t*sqrt(norm(gd1,2)^2 + norm(gd2,2)^2)
        t = beta*t;
        new_angs = angs - t*gd1;
        new_Bs = Bs - t*gd2;
        new_res = residual(node,u,omega,x0,y0,c0,new_angs,new_Bs);
    end
    
    t = t
    new_res = new_res
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
        
        rd = 40*h;
        M = 32;
        angl = linspace(0,2*pi,M+1) ; 
        ang = angl(1:M) ;
        x = x0 + rd*cos(ang');
        y = y0 + rd*sin(ang');
        [~, x, y] = xy_to_index(node,x,y);
        
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
        
        rd = 40*h;
        M = 32;
        angl = linspace(0,2*pi,M+1) ; 
        ang = angl(1:M) ;
        x = x0 + rd*cos(ang');
        y = y0 + rd*sin(ang');
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


