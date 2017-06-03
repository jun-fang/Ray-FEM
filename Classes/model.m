classdef model < handle
    % The model class is a encapsulated model with all the paramaters
    % needed to solve the Helmholtz equation inside it
    
    properties
        node  % the node structure
        elem 
        omega 
        npml    % number of points of the PML
        wpml
        sigmaMax
        speed
        fquadorder
        local_solver % local solver encapsulated in here
        plt
        H       % Helmholtz operator
        freeNode
        %bdyConditions
    end
    
    methods
        % function to build the model and save the data needed to solve it
        function init(obj, node,elem,omega,wpml,sigmaMax,speed,fquadorder)
            
            % assembling the full matrix
            A = assemble_Helmholtz_matrix(node, elem, omega, wpml, ...
                                          sigmaMax, speed, fquadorder);
            
            obj.node = node;
            obj.elem = elem; 
            obj.omega = omega;
            %obj.npml    % number of points of the PML
            obj.wpml = wpml;
            obj.sigmaMax = sigmaMax;
            obj.speed = speed ;
            obj.fquadorder = fquadorder;
            % Boundary conditions for the PML Case
            [~,~,isBdNode] = findboundary(elem);
            obj.freeNode = find(~isBdNode);
            % we need to add all the local solver properties
            % in particulat which are the indices of the degrees of 
            % freedom at the boundary
            
            % Imposing the homogeneous Dirichlet boundary conditions
            fprintf('Assembling the Helmholtz matrix \n');
            obj.H = A(obj.freeNode,obj.freeNode);
            
            %obj.local_solver = solver(obj.H);
        end
         
        function initRay(obj, node,elem,omega,wpml,sigmaMax,speed,ray,fquadorder)
            % init system using the ray based matrices
            % assembling the full matrix
            A = assemble_Helmholtz_matrix_RayFEM(node, elem, omega, wpml, ...
                                          sigmaMax, speed, ray, fquadorder);
            
            obj.node = node;
            obj.elem = elem; 
            obj.omega = omega;
            %obj.npml    % number of points of the PML
            obj.wpml = wpml;
            obj.sigmaMax = sigmaMax;
            obj.speed = speed ;
            obj.fquadorder = fquadorder;
            % Boundary conditions for the PML Case
            [~,~,isBdNode] = findboundary(elem);
            obj.freeNode = find(~isBdNode);
            % we need to add all the local solver properties
            % in particulat which are the indices of the degrees of 
            % freedom at the boundary
            
            % Imposing the homogeneous Dirichlet boundary conditions
            fprintf('Assembling the Helmholtz matrix \n');
            obj.H = A(obj.freeNode,obj.freeNode);
            
            %obj.local_solver = solver(obj.H);
        end
        % TODO extend the definition of the function in order to handle 
        % boundary conditions. 
        function u = solve(obj, f)
        %  u = solve(obj, f)
        % function to solve the Helmholtz equation
        % input:  f source, this can be a either a vector or a function
        %                   handle
           if  isa(f,'function_handle')
               u = zeros(size(obj.node,1),1);
               % f is a function handle, then we need to project it into
               % the Galerking space
               b = assemble_RHS(obj.node,obj.elem,f,obj.fquadorder);
               % we use only the interior points 
               % and we call the solve function in the solver class
               u(obj.freeNode) = obj.local_solver.solve(b(obj.freeNode));
           elseif isvector(f)
               % if f is a vector
               u = zeros(size(obj.node,1),1);
               u(obj.freeNode) = obj.local_solver.solve(f(obj.freeNode));
           end
        end
        
        function u = solveInt(obj, f)
        %  u = solveInt(obj, f)
        % function to solve the Helmholtz equation for only the interior 
        % points, this is necessary for the fast solver (otherwise is a
        % huge pain 
        % input:  f source, this can be a either a vector or a function
        %                   handle
           if isvector(f)
               % if f is a vector
               u = obj.local_solver.solve(f);
           end
        end
        
        function LU_factorization(obj)
            % defining the local factorization, which is encapsulated
            % inside the local_solver (see the 
            obj.local_solver = solver(obj.H); 
        end
        
        
       function showresult(obj,u,viewangle)
        %% SHOWRESULT display the mesh and the solution 
        %
        %  showresult(node,elem,u,viewangle) displays the mesh and the solution in
        %  one figure. The left one is the mesh, the middle
        %  one is the contour of the solution, and the right one is the graph of
        %  the function. The last viewangle is used to adjust the view angle of the
        %  graph of the function.
        %

        if (length(u) == size(obj.elem,1)) || (length(u) == size(obj.node,1))
            % show mesh
            %set(gcf,'Units','normal'); 
            %set(gcf,'Position',[0.25,0.25,0.45,0.25]); 
            subplot(1,2,1); 
            showsolution(obj.node,obj.elem,u,2);
            axis 'equal'
            colorbar;
            %     pause(0.05)
            subplot(1,2,2); 
            showsolution(obj.node,obj.elem,u);
            axis 'equal'
            if nargin>2
                view(viewangle);
            end
            %     pause(0.05)
        else
        showmesh(obj.node,obj.elem);
        end
       end
        
    end  
    
end

