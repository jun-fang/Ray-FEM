classdef subdomain < handle
    % The model class is a encapsulated model with all the paramaters
    % needed to solve the Helmholtz equation inside it
    
    properties
        model
        indIntLocal
        indIntGlobal  
        indxn
        indxnp
        indx0
        indx1
        n1   % number of degreed of freedom at the top interface
        nn   % number of degrees of freedom at the bottom interface
        position
        %bdyConditions
    end
    
    methods
        %TODO find the correct data for the subdomain
        function initModel(obj, node,elem,omega,wpml,sigmaMax,pde,fquadorder)
            
            % assembling the full matrix
            A = assemble_Helmholtz_matrix(node, elem, omega, wpml, ...
                                          sigmaMax, pde.speed, fquadorder);
            
            obj.node = node;
            obj.elem = elem; 
            obj.omega = omega;
            %obj.npml    % number of points of the PML
            obj.wpml = wpml;
            obj.sigmaMax = sigmaMax;
            obj.pde = pde ;
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
            
            obj.local_solver = solver(obj.H);
        end
        
        function init_Model(obj,model,indIntLocal, indIntGlobal, indxn, indxnp, indx0, indx1, position)
            % function to initialize the subdomain given that hte model has
            % already been computed.
            obj.model = model;
            obj.indIntLocal = indIntLocal;
            obj.indIntGlobal  = indIntGlobal;
            obj.nn = length(indxn);
            obj.indxn  = indxn;
            obj.indxnp = indxnp;
            obj.n1 = length(indx0);
            obj.indx0  = indx0;
            obj.indx1  = indx1;
            obj.position = position;
        end
        
        function u = solveTraces(obj, f, v0,v1, vn,vnp)
        % u = solveInt(obj, f)
        % function to solve the Helmholtz equation for only the interior 
        % points, this is necessary for the fast solver (otherwise is a
        % huge pain 
        % input:    f vector encoding the source
        %           v0 trace of v at 0
        %           v1 trace of v at 1
        %           vn trace of v at n
        %           vnp trace of v at np
        
        % size checks
        
        if obj.position ~= 'N'
            if (length(v0) ~= obj.n1 ||  length(v1) ~= obj.n1)
                error('vector size of the traces 0 or 1 are not consistent')
            end
            
            f(obj.indx0,:) = f(obj.indx0,:) + obj.model.H(obj.indx0,obj.indx1)*v1;
            f(obj.indx1,:) = f(obj.indx1,:) - obj.model.H(obj.indx1,obj.indx0)*v0;
        end
        
        
        if obj.position ~= 'S'
        if (length(vn) ~= obj.nn ||  length(vnp) ~= obj.nn)
            error('vector size of the traces 0 or 1 are not consistent')
        end
        
            f(obj.indxn,:) = f(obj.indxn,:)  - obj.model.H(obj.indxn ,obj.indxnp)*vnp;
            f(obj.indxnp,:)= f(obj.indxnp,:) + obj.model.H(obj.indxnp,obj.indxn )*vn;
        
        end
        
        u = obj.model.solveInt(f);
        
        end
        
        function [u0, u1, un, unp] = applyBlockOperator(obj, v0,v1, vn,vnp)
            % function to apply the integral operator to the traces of the
            % wavefield at the interfaces.
            f = zeros(size(obj.model.node(obj.model.freeNode,1),1),1);
            u = obj.solveTraces(f, v0,v1, vn,vnp);
            [u0, u1, un, unp] = obj.extractLocalTraces(u);
        end
        
        function [u0, u1, un, unp] = extractLocalTraces(obj,u)
            % TODO put sizes checks everywhere!!! 
            u0  = u(obj.indx0);
            u1  = u(obj.indx1);
            un  = u(obj.indxn);
            unp = u(obj.indxnp);
        end
        
        function LU_factorization(obj)
            % defining the local factorization, which is encapsulated
            % inside the local_solver (see the 
            obj.model.LU_factorization(); 
        end
        
        
        function fIntLocal = extractLocalSource(obj, f)
            % function to extract the information from the global source to
            % to a local one.
            % create a local vector with the correct sizes (carefull with
            % the extra degrees of freedom from the boudnary that have to
            % be removes
            fIntLocal = zeros(size(obj.model.node(obj.model.freeNode,1),1),1);
            % copying from the global source to the local one
            fIntLocal(obj.indIntLocal) = f(obj.indIntGlobal);
        end
        
        function uIntLocal = extractLocalInt(obj, u)
            % function to extract the wavefield information from the
            % interios degrees of freedom of a local wavefield
            uIntLocal = u(obj.indIntLocal); 
        end
        
        % TODO extend the definition of the function in order to handle 
        % boundary conditions. 
        function u = solveLocal(obj, f)
        %  u = solve(obj, f)
        % function to solve the Helmholtz equation
        % input:  f source, this can be a either a vector or a function
        %                   handle
         u = obj.model.solve(f);
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
               u = obj.model.solveInt(f);
           end
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
            colorbar;
            %     pause(0.05)
            subplot(1,2,2); 
            showsolution(obj.node,obj.elem,u);
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

