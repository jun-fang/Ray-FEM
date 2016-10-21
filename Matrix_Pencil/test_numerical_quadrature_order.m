function test_numerical_quadrature_order(omega,h,a,b)

% a reference triangle
node = [0 0; 1 0; 0 1]*h;
elem = [1 2 3];
area = 0.5*h*h;
err = zeros(9,2);
A = 1i*omega*a;
B = 1i*omega*b;
for n = 1:9    
    % get quadrature points
    [lambda,weight] = quadpts(n);
    nQuad = size(lambda,1);
    t1 = 0; 
    for p = 1:nQuad
        % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        t1 = t1 + weight(p)*f1(pxy(1),pxy(2),omega,a,b);
    end                
    t1 = t1*area;
    err(n,1) = omega^2*abs(t1 - f2(h,A,B));
end
display('Table: Error of quadrature for e^{i \omega d \cdot x}')
colname = {'n','Error'};
disptable(colname,(1:9)',[],err(:,1),'%0.5e');
end

function z = f1(x,y,omega,a,b)
z = exp(1i*omega*(a*x + b*y));
end

function z = f2(h,A,B)
z = (-A*exp(B*h) + exp(A*h)*B + A-B)/(A^2*B-A*B^2);
end

%% Results
% Let T be the triangle formed by (0,0), (1,0), and (0,1). 
%
% Error1 is for the integral 
% $\int _{T} x^n + y^n \, dxdy$. 
% It should be numerically exact. 
%
% Error2 is for the integral 
% $\int _{T} \sin(x+y) \, dxdy$.
% It decays as n increas.
%
% See the doc for qudrature rules in <matlab:ifemdoc('quadpts') quadpts>.
