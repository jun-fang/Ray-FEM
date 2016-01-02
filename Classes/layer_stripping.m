classdef layer_stripping < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    % what are the properties that we need? 
    
    
    properties
        M           % in this structure we are going to save all the models
        M_base      % this is the total basic model for comparaison / must be deleted in a newer version
        overlap     % overlap between all the models
        n_layer     % number of layers
        npml        % number of points for the pml's
        c           % the full velocity model
        order       % order of accuracy
        h           % the discretization step
        sigmaMax    % maximum value for the PML profile
        K           % wavevectir for the helmholtz equation
        omega       % frequency for the helmholt equation
        y_sample    % we need to know the sampling depth
        x_extrp     % the extrapolation depth
        nxi     % number of interior points in x
        nyi     % number of interior points in y
        ni      % total number of poitn in the interior grid
        xi      %interior grid in x
        yi      %interior grid in x
        Xi      %interior 2_D grid in x
        Yi      %interior 2_D grid in y
        nx  
        ny 
        x
        y
        X      % 2D grid of the whole model
        Y      % 2D grid of the whole model
        n      % total number of interior points
    end % properties
    
    methods
        function init(obj , c, n_layer, npml, overlap, K, sigmaMax, h, order)
            obj.c = c;
            obj.n_layer = n_layer; 
            obj.npml    = npml; 
            obj.overlap = overlap; 
            obj.K       = K;
            obj.omega   = 2*pi*K;
            obj.sigmaMax= sigmaMax;
            obj.h       = h;
            obj.order   = order;
            obj.nxi     = size(c,2);     
            obj.nyi     = size(c,1);     
            obj.ni      = obj.nxi*obj.nyi;
            obj.nx      = obj.nxi + 2*npml;
            obj.ny      = obj.nyi + 2*npml;
            obj.n       = obj.nx*obj.ny;
            
            local_ni = floor(obj.nyi/obj.n_layer); 
            
            obj.M   = {};
%             obj.M_base = model; 
%             obj.M_base.init(c,h, npml, K, sigmaMax, order, [0 0]); % we suppose that we always start from zero elevation
            
            obj.xi = h*(0:obj.nxi-1);
            obj.yi = h*(0:obj.nyi-1);
            obj.x  = [obj.xi(1)+(-npml:-1)*h obj.xi obj.xi(end)+(1:npml)*h];
            obj.y  = [obj.yi(1)+(-npml:-1)*h obj.yi obj.yi(end)+(1:npml)*h]; 

            [obj.Xi,obj.Yi] = meshgrid(obj.xi,obj.yi);
            [obj.X ,obj.Y ] = meshgrid(obj.x ,obj.y );
            
            % although it's called x_array, the array is in depth not in
            % offset
            ind_array  = 1:local_ni:obj.nyi;
            ind_array  = [ ind_array(1:obj.n_layer) obj.nyi+1];
            %ind_array(end) = obj.nyi+1; % because we count from zero
            
            if n_layer > 2 % setting the "overlap" to have different boundaries 
                ind_array(2:end-1) = ind_array(2:end-1)+ overlap;
            end
            %the overlap is used to have different set of boundaries
            
            %this are 1D array and they are orderred
            if n_layer > 1
                y_sampling_down = obj.yi(ind_array(1:end-1));
                x_extrp_down    = obj.yi(ind_array(2:end)-1);
            else
                y_sampling_down = obj.yi(1);
                x_extrp_down    = obj.yi(end);
            end
            obj.y_sample = y_sampling_down;
            obj.x_extrp  = x_extrp_down;
            
            c_array = {}; 
            c = reshape(c,obj.nyi,obj.nxi); 
            for ii = 1:obj.n_layer
                c_aux = c(find( (obj.Yi>= (y_sampling_down(ii)-10*eps)).*(obj.Yi<=(x_extrp_down(ii)+10*eps))) );
                % have to reshape c, cuz the model init uses the geomtry of
                % c to initialize the grid
                c_aux = reshape(c_aux, numel(c_aux)/obj.nxi, obj.nxi);
                c_array{ii} = c_aux;
            end
            
            for ii = 1:1:obj.n_layer                
                obj.M{ii} = model; 
                obj.M{ii}.init( c_array{ii} ,h, npml, K, sigmaMax, order, [0 obj.y_sample(ii)] );
                % have to specify a flag here in order to show this
                %figure(30 + ii); clf;
                %DisplayField( obj.M{ii}.c, obj.M{ii}.xi, obj.M{ii}.yi); 
                % gotta be carefull of how this is defined 
            end
        end
        
        function factorization_slabs(obj)
            for ii = 1:1:obj.n_layer                           
                obj.M{ii}.LU_factorization();
            end
        end
        
        % function to copy information from one object to another
        % ideal for debugging, we only need to initialize one object save
        % it and then copy it to the new one :)       
        function copy(obj, obj_2)
            % have to put safe guard that we are in fact copying layered
            % media objects and not something else.
            obj.M           = obj_2.M;
            obj.M_base      = onj_2.M_base;
            obj.overlap     = obj_2.overlap;
            obj.n_layer     = obj_2.n_layer;
            obj.npml        = obj_2.npml;
            obj.c           = obj_2.c ;
            obj.order       = obj_2.order;
            obj.h           = obj_2.h  ;
            obj.sigmaMax    = obj_2.sigmaMax ;
            obj.K           = obj_2.K;
            obj.omega       = obj_2.omega    ;
            obj.y_sample    = obj_2.y_sample  ;
            obj.x_extrp     = obj_2.x_extrp ;
            obj.nxi         = obj_2.nxi  ;
            obj.nyi         = obj_2.nyi   ;
            obj.ni          = obj_2.ni  ;
            obj.xi          = obj_2.xi ;
            obj.yi          = obj_2.yi ;
            obj.Xi          = obj_2.Xi ;
            obj.Yi          = obj_2.Yi  ;
            obj.nx          = obj_2.nx ;
            obj.ny          = obj_2.ny;
            obj.x           = obj_2.x;
            obj.y           = obj_2.y;
            obj.X           = obj_2.X ;
            obj.Y           = obj_2.Y  ;
            obj.n           = obj_2.n   ;
        end
        
        % this fucntion will extract the traces of the solution at the
        % ii-th interface
        function [u_trace, Dx_u_trace]= extract_interface_trace(obj,u, ii)   
             u_trace = u(find(obj.Y== obj.x_extrp(ii)));      
             Dx_u_trace = (u(find(abs((obj.Y-(obj.x_extrp(ii)+obj.h)))< 10*eps)) - ... 
                           u(find(abs((obj.Y-(obj.x_extrp(ii)-obj.h)))< 10*eps)) )/(2*obj.h) ;   
        end
        
        function [u_trace, u_trace_p]= extract_local_dirichlet_trace(obj,u, ii, arg)
            % by default it computes the derivative at the extrapolation
            % line, arg is a cell of strings
         
            % downgoing sweep with a descentered stencil of first order
            if strcmp(arg, 'down')
                u_trace    = u(find(abs(obj.M{ii}.Y- (obj.x_extrp(ii)))<10*eps));
                u_trace_p = u(find(abs((obj.M{ii}.Y-( obj.x_extrp(ii) + obj.h)))< 10*eps)) ;
            % upgoing sweep with a descentered stencil of first order
            elseif strcmp(arg, 'up')
                u_trace    = u(find(abs(obj.M{ii}.Y- (obj.y_sample(ii)))<10*eps));
                u_trace_p = u(find(abs((obj.M{ii}.Y-(obj.y_sample(ii) - obj.h)))< 10*eps)) ;
            end           
        end 
        
        
        
        function u_trace = extract_trace(obj,u, y_depth)   
             u_trace = u(find(abs(obj.Y - y_depth) < 2*eps ));      
        end
        
        % this function extracts the trace and the derivative of the local
        % solution 
        function [u_trace, Dx_u_trace]= extract_interface_local_trace(obj,u, ii, arg)
            % by default it computes the derivative at the extrapolation
            % line, arg is a cell of strings
         
            % downgoing sweep with a descentered stencil of first order
            if strcmp(arg, 'down')
                u_trace = u(find(abs(obj.M{ii}.Y- (obj.x_extrp(ii)))<10*eps));
                Dx_u_trace = ( - u_trace + ...
                    u(find(abs((obj.M{ii}.Y-( obj.x_extrp(ii) + obj.h)))< 10*eps)))/(obj.h) ;
            % upgoing sweep with a descentered stencil of first order
            elseif strcmp(arg, 'up')
                u_trace = u(find(abs(obj.M{ii}.Y- (obj.y_sample(ii)))<10*eps));
                Dx_u_trace = (  u_trace - ...
                    u(find(abs((obj.M{ii}.Y-(obj.y_sample(ii) -obj.h)))< 10*eps)))/(obj.h) ;
            end           
        end 
        
        
        
        % this function is made to split the sources in such a way that it
        % can be handled by every layer. 
        function s_array =  source_splitting(obj,s, flag)
            if nargin == 2
                s_array = {};
                for ii = 1:1: obj.n_layer
                    s_aux = s(find( (obj.Y>=obj.y_sample(ii)).*(obj.Y<=obj.x_extrp(ii))) );
                    s_aux = reshape(s_aux, numel(s_aux)/obj.nx, obj.nx);
                    % padding zeros up and down (this can be a problem if the
                    % source is no completely inside the domain
                    s_array{ii} = [ zeros(obj.npml,obj.nx); s_aux ; zeros(obj.npml,obj.nx) ];
                end
            elseif nargin == 3
                if strcmp(flag,'extended') % for the single partition test
                    s_array = {};
                    for ii = 1:1: obj.n_layer
                        % we have to split the source and add one more pixel up
                        % and down
                        s_aux = s(find( (obj.Y>= obj.y_sample(ii)-obj.h-10*eps).*(obj.Y<=obj.x_extrp(ii)+obj.h+10*eps)));
                        s_aux = reshape(s_aux, numel(s_aux)/obj.nx, obj.nx);
                        % padding zeros up and down (this can be a problem if the
                        % source is no completely inside the domain
                        s_array{ii} = [ zeros(obj.npml-1,obj.nx); s_aux ; zeros(obj.npml-1,obj.nx) ];
                    end
                end
            end
                
            end % function source_splitting
            
        
        % Function that perform the local solves of every sub-domain
        function u_array =  solve_local(obj,s, flag)
            % have to implement the solver with Dirichelt boundary conditions as well
            
            u_array = {};
            if nargin == 2
                s_array = source_splitting(obj, s);
            else
                s_array = source_splitting(obj, s, flag); % for the single partition test
            end
            
            parfor ii = 1:obj.n_layer
                u_aux       = obj.M{ii}.solve(s_array{ii});
                u_aux       = reshape(u_aux, obj.M{ii}.ny, obj.M{ii}.nx);
                u_array{ii} = u_aux;
            end 
                
        end % function solve_local
        
        % Function to concatenate a series of partial local solution in
        % order to obtain a global solution. Its main application is to
        % compare different solutions 
        
        function u_global = concatenate(obj, u_array, flag)
            
            aux_counter = 0;
            if nargin==2
                u_global = zeros(size(obj.Xi));
                for ii = 1:numel(u_array)
                    % we should have continuity in the overlaping interfacing
                    % so this should be more than enough. Otherwise I need to
                    % change this in order to consider the jump conditions.
                    
                    u_reordered = reshape(u_array{ii}, obj.M{ii}.ny, obj.M{ii}.nx);
                    u_global( aux_counter+1:obj.M{ii}.nyi+aux_counter, 1:obj.M{ii}.nxi) = ...
                        u_reordered( obj.M{ii}.npml+1:obj.M{ii}.npml+obj.M{ii}.nyi, obj.M{ii}.npml+1:obj.M{ii}.npml+obj.M{ii}.nxi );
                    aux_counter = aux_counter + (obj.M{ii}.nyi) ;
                end
                
            elseif nargin==3 && strcmp(flag,'extended') % this flag is to add the PML 
                u_global = zeros(size(obj.X));
                for ii = 1:numel(u_array)
                    % we should have continuity in the loverlaping interfacing
                    % so this should be more than enough. Otherwise I need to
                    % change this in order to consider the jump conditions.
                    
                    % have to reorder the solutions (this is legacy from
                    % Vincent;s code)
                    u_reordered = reshape(u_array{ii}, obj.M{ii}.ny, obj.M{ii}.nx);
                    if ii~= 1 && ii~= numel(u_array)
                        
                        u_global( aux_counter+1:obj.M{ii}.nyi+aux_counter, 1:obj.M{ii}.nx) = ...
                            u_reordered( obj.M{ii}.npml+1:obj.M{ii}.npml+obj.M{ii}.nyi, 1:obj.M{ii}.nx);
                        aux_counter = aux_counter + obj.M{ii}.nyi  ;
                        
                    elseif ii == 1 % in the first one we need to add the pml at the begining
                        
                        u_global( aux_counter+1:obj.M{ii}.nyi+aux_counter+obj.M{ii}.npml, 1:obj.M{ii}.nx) = ...
                            u_reordered(1:obj.M{ii}.npml+obj.M{ii}.nyi, 1:obj.M{ii}.nx);
                        aux_counter = aux_counter + obj.M{ii}.npml + obj.M{ii}.nyi ;
                        
                    elseif ii == numel(u_array) % last one, have to add the pml as well
                        
                        u_global( aux_counter+1:obj.M{ii}.nyi+aux_counter+obj.M{ii}.npml, 1:obj.M{ii}.nx) = ...
                            u_reordered( obj.M{ii}.npml+1:2*obj.M{ii}.npml+obj.M{ii}.nyi, 1:obj.M{ii}.nx);
                        aux_counter = aux_counter + obj.M{ii}.nyi ;
                        
                    end
                end
                
                
            end
            
            
            
        end % fucntion concatenate
        
        % we compute an array of green's function objects to compute the
        % extrapolation. all the information is already contained in the
        % layered media class
        % need to define upwards and downwatrds GReen's fucntions
        
         function GG = computing_Green_functions_optimized(obj, arg)
            %by default the sweep is downwards if we put arg then it will
            %differentiate beetween up and down
            if nargin==2
                % this is the factorization going upwards have to be
                % carefull on the signs and the normal is defined
                if strcmp(arg{1}, 'up')
                    GG = {};
                    for ii = 1: obj.n_layer
                        GG{ii} = Green_functions;
                    end
                    if numel(arg) == 1
                        for ii = 1: obj.n_layer
                            fprintf('Extracting the upwards Green function for the %i th layer \n', ii);
                            GG{ii}.init(obj.M{ii}.local_solver, obj.M{ii}.X, obj.M{ii}.Y, [obj.x(1) obj.x(end)], ...
                                [obj.x_extrp(ii) obj.y_sample(ii)], -1);
                        end
                    elseif numel(arg) == 3
                        for ii = 1: obj.n_layer
                            fprintf('Extracting the upwards Green function for the %i th layer and compressing them as PLR matrices \n', ii);
                            GG{ii}.init(obj.M{ii}.local_solver, obj.M{ii}.X, obj.M{ii}.Y, [obj.x(1) obj.x(end)], ...
                                [obj.x_extrp(ii) obj.y_sample(ii)], -1, [arg{2} arg{3}] );
                        end
                    end
                elseif   strcmp(arg{1}, 'down')
                    GG = {};
                    for ii = 1: obj.n_layer
                        GG{ii} = Green_functions;
                    end
                    if numel(arg) == 1
                        for ii = 1: obj.n_layer
                            fprintf('Extracting the downwards Green function for the %i th layer \n', ii);
                            GG{ii}.init(obj.M{ii}.local_solver, obj.M{ii}.X, obj.M{ii}.Y, [obj.x(1) obj.x(end)], ...
                                [obj.y_sample(ii) obj.x_extrp(ii)], 1 );
                        end
                    elseif numel(arg) == 3 % for compressed matrices
                        for ii = 1: obj.n_layer
                            fprintf('Extracting the downwards Green function for the %i th layer and compressing them as PLR matrices \n', ii);
                            GG{ii}.init(obj.M{ii}.local_solver, obj.M{ii}.X, obj.M{ii}.Y, [obj.x(1) obj.x(end)], ...
                                [obj.y_sample(ii) obj.x_extrp(ii)],1, [arg{2} arg{3}] );
                        end
                    end
                end
                
                
            end
            
        end % function computing_Green_function
        
        % function to display the local wavefields given by solving the
        % local Systems
        function DisplayField(obj, u_array)
            % have to check the lenght of an array 
            for ii = 1:obj.n_layer
                figure(100+ii);
                DisplayField(u_array{ii}, obj.M{ii}.x, obj.M{ii}.y, obj.M{ii}.npml);
            end
        end % function DisplayField
        
        %small fucntion that converts the sparse data on a 
        function s = struct2field(obj, s_struct)
            s = zeros(size(obj.X));
            for ii = 1:numel(s_struct.y)
               ind_y = find(abs(obj.Y-s_struct.y(ii))<2*eps);
               s(ind_y) = s_struct.u(:,ii);
            end
            
        end
            
    end % methods
    
end % classdef

 