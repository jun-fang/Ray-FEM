classdef layer_stripping < handle
    % layer_stripping: class that contains the information of the whole
    % solver. It contains the data structures fo the horizontal layers and
    % which contain the information of each subdomain 
    
    properties
        M           % in this structure we are going to save all the models
        MBase      % this is the total basic model for comparaison / must be deleted in a newer version
        nLayer     % number of layers
        P
    end % properties
    
    methods
        
        function init_model(obj, M, Mbase, nLayer)
            obj.M = M; 
            obj.MBase= Mbase; 
            obj.nLayer = nLayer;
        end
        
        % factorizing the Helmholtz matrices at within each slab 
        function factorization_slabs(obj)
            for ii = 1:1:obj.nLayer                           
                obj.M{ii}.model.LU_factorization();
            end
        end % function factorization_slabs
                 
        function [uy0,uy1, uyN, uyNp ]= extractLocalTrace(obj,u, ii)
            % this function extracts the boundary traces of u, given that u
            % is the local solution at the ii-th subdomain
            
             [uy0,uy1, uyN, uyNp ] = obj.M{ii}.extractLocalTraces(u);
                  
        end % function extractLocalTrace
        
        function uGamma = extractGlobalTraces(obj, u)
        % function to extract the traces from a global solution
            uSplit =  obj.source_splitting(u);
            
           uGamma = [];
           ii = 1;
           [~,~, uyN, ~ ]= obj.extractLocalTrace(uSplit{ii}, ii);
           uGamma = [uGamma; uyN];
           
           for ii = 2:obj.nLayer-1
                [~,uy1, uyN, ~ ]= obj.extractLocalTrace(uSplit{ii}, ii);
                uGamma = [uGamma; uy1; uyN];
           end
           ii = obj.nLayer;
           [~,uy1, ~, ~ ]= obj.extractLocalTrace(uSplit{ii}, ii);
           uGamma = [uGamma; uy1];    
        
        end
                
        % Function that perform the local solves of every sub-domain
        function u_array =  solve_local(obj,s, flag, flagpol)
            % have to implement the solver with Dirichelt boundary conditions as well
            
            u_array = {};
            if nargin == 2 || ~strcmp(flag,'extended')  
                s_array = obj.source_splitting( s);
            else 
                s_array = obj.source_splitting( s, flag); % for the single partition test
            end
            
            
            for ii = 1:obj.nLayer
                if nargin == 4
                    u_aux       = obj.M{ii}.solveInt(s_array{ii}, flagpol);
                else
                    u_aux       = obj.M{ii}.solveInt(s_array{ii});
                end
                u_array{ii} = u_aux;
            end
            
        end % function solve_local
        
        function f = formRHS(obj,uArray)
           % function to form the RHS for the basic formulation (no
           % polarization)
           ii = 1;
           f = [];
           [~,~, uyN, ~ ]= extractLocalTrace(obj,uArray{ii}, ii);
           f = [f; uyN];
           
           for ii = 2:obj.nLayer-1
                [~,uy1, uyN, ~ ]= extractLocalTrace(obj,uArray{ii}, ii);
                f = [f; uy1; uyN];
           end
           ii = obj.nLayer;
           [~,uy1, ~, ~ ]= extractLocalTrace(obj,uArray{ii}, ii);
           f = [f; uy1];    
        end % function formRHS
        
        function ff = formPolarizedRHS(obj,uArray)
           % function to form the RHS for the basic formulation (no
           % polarization)
           ii = 1;
           f = [];
           f0 = [];
           [~,~, uyN, uyNp ]= extractLocalTrace(obj,uArray{ii}, ii);
           f = [f; uyN];
           f0 = [f0; uyNp];
           
           for ii = 2:obj.nLayer-1
                [uy0,uy1, uyN, uyNp]= extractLocalTrace(obj,uArray{ii}, ii);
                f = [f; uy1; uyN];
                f0 = [f0; uy0; uyNp];
           end
           ii = obj.nLayer;
           [uy0,uy1, ~, ~ ]= extractLocalTrace(obj,uArray{ii}, ii);
           f = [f; uy1];  
           f0 = [f0; uy0];
           
           ff = [f; f0];
           
        end % function formPolarizedRHS
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % functions that apply the integral system (normal or polarized)
        % and the different decompositions in matrix free form, relying on
        % the local solve to accomplish this.
        
        function Mu = applyM(obj, uBdy)
            % function that applies M, the integral matrix in matrix free form, it
            % only relies on local solves to apply the matrix.
            
            Mu = zeros(size(uBdy));
            
            ii = 1;
            % parsing the information at the boundary
            indlocal1 = (1:obj.M{ii}.n1);  % defining the local indeces
            indlocaln = (1:obj.M{ii}.nn);
            indStart = 0; % defining the last index
            indShift1 = obj.M{ii}.n1;
            indShiftn = obj.M{ii}.nn;
            % parsing the vector uBdy
            uy0  = 0*uBdy(indlocal1);
            uy1  = 0*uBdy(indlocal1);
            uyN  = uBdy(          indlocaln);
            uyNp = uBdy(indShiftn+indlocaln);

            [vy0, vy1, vyN, vyNp] = obj.M{1}.applyBlockOperator( uy0, uy1, uyN, uyNp );
            Mu(indlocaln) = vyN - uyN;
            
            for ii = 2:obj.nLayer-1
                % parsing the information at the boundary
                % we need to remember that in the general case, we won't
                % have the same amount of interface data at each layer
                % this part should be modified accordingly
                indlocal1 = indStart + (1:obj.M{ii}.n1);  % defining the local indeces
                indlocaln = indStart + 2*obj.M{ii}.n1 + (1:obj.M{ii}.nn);
                indStart =  max(indlocal1) + obj.M{ii}.n1;
                indShift1 = obj.M{ii}.n1;
                indShiftn = obj.M{ii}.nn;
                
               
                % parsing the vector uBdy
                uy0  = uBdy(          indlocal1);
                uy1  = uBdy(indShift1+indlocal1);
                uyN  = uBdy(          indlocaln);
                uyNp = uBdy(indShiftn+indlocaln);
                
                [vy0, vy1, vyN, vyNp]   = obj.M{ii}.applyBlockOperator( uy0, uy1, uyN, uyNp );
                Mu(indShift1+indlocal1) = vy1 - uy1;
                Mu(          indlocaln) = vyN - uyN;
                
            end
            
            ii = obj.nLayer;
            indlocal1 = indStart + (1:obj.M{ii}.n1);  % defining the local indeces
            indlocaln = indStart + 2*obj.M{ii}.n1 + (1:obj.M{ii}.nn);
            indStart =  max(indlocal1) + obj.M{ii}.n1;
            indShift1 = obj.M{ii}.n1;
            indShiftn = obj.M{ii}.nn;
            % parsing the vector
            uy0  = uBdy(         indlocal1);
            uy1  = uBdy(indShift1+indlocal1);
            uyN  = 0*uBdy(indlocal1);
            uyNp = 0*uBdy(indlocal1);
            
            [vy0, vy1, vyN, vyNp] = obj.M{ii}.applyBlockOperator( uy0, uy1, uyN, uyNp );
            Mu(indShift1+indlocal1)          = vy1 - uy1;
            
        end % function applyM
        
        function Mu = applyMM1(obj,uBdy)
            % function optimized to compute the application of Mup and
            % Mdown simultaneously, 
            uDown = uBdy(1:end/2);
            uUp   = uBdy(end/2+1:end);
            
            Mu = zeros(length(uBdy)/2,1);
            
            ii = 1;
            % parsing the information at the boundary
            indlocal1 = (1:obj.M{ii}.n1);  % defining the local indeces
            indlocaln = (1:obj.M{ii}.nn);
            indStart = 0; % defining the last index
            indShift1 = obj.M{ii}.n1;
            indShiftn = obj.M{ii}.nn;
            % parsing the vector uBdy
            uy0Up  = 0*uUp(indlocal1);
            uy1Up  = 0*uUp(indlocal1);
            uyNUp  = uUp(          indlocaln);
            uyNpUp = uUp(indShiftn+indlocaln);
            uyNDown  = uDown(          indlocaln);
            uyNpDown = uDown(indShiftn+indlocaln);

            [~, ~, vyN, ~] = obj.M{1}.applyBlockOperator(uy0Up,...
                                             uy1Up, uyNUp, uyNpUp);
            Mu(indlocaln) = vyN - uyNUp - uyNDown;
            
            for ii = 2:obj.nLayer-1
                % parsing the information at the boundary
                indlocal1 = indStart + (1:obj.M{ii}.n1);  % defining the local indeces
                indlocaln = indStart + 2*obj.M{ii}.n1 + (1:obj.M{ii}.nn);
                indStart =  max(indlocal1) + obj.M{ii}.n1;
                indShift1 = obj.M{ii}.n1;
                indShiftn = obj.M{ii}.nn;
                
               
                 % parsing the vector uBdy
                uy0Down  = uDown(          indlocal1);
                uy1Down  = uDown(indShift1+indlocal1);
                uyNDown  = uDown(          indlocaln);
                uyNpDown = uDown(indShiftn+indlocaln);
                
                  % parsing the vector uBdy
                uy0Up  = uUp(          indlocal1);
                uy1Up  = uUp(indShift1+indlocal1);
                uyNUp = uUp(          indlocaln);
                uyNpUp = uUp(indShiftn+indlocaln);
                
                                        
                [~, vy1, ~, ~] = obj.M{ii}.applyBlockOperator( uy0Down, uy1Down, uyNDown+uyNUp, uyNpDown +uyNpUp );
                Mu(indShift1+indlocal1) = vy1 - uy1Down - uy1Up;
                
                [~, ~, vyN, ~]   = obj.M{ii}.applyBlockOperator( uy0Up+uy0Down, uy1Up+uy1Down, uyNUp, uyNpUp );
                Mu(          indlocaln) = vyN - uyNUp - uyNDown;
                
                
            end
            
            ii = obj.nLayer;
            indlocal1 = indStart + (1:obj.M{ii}.n1);  % defining the local indeces
            indlocaln = indStart + 2*obj.M{ii}.n1 + (1:obj.M{ii}.nn);
            indStart =  max(indlocal1) + obj.M{ii}.n1;
            indShift1 = obj.M{ii}.n1;
            indShiftn = obj.M{ii}.nn;
            % parsing the vector
            uy0Down  = uDown(         indlocal1);
            uy1Down  = uDown(indShift1+indlocal1);
            uyNDown  = 0*uDown(indlocal1);
            uyNpDown = 0*uDown(indlocal1);
            
            [~, vy1, ~, ~] = obj.M{ii}.applyBlockOperator( uy0Down, uy1Down, 0*uyNDown, 0*uyNpDown );
            Mu(indShift1+indlocal1)          = vy1 - uy1Down - uy1Up ;
            
            
        end
        
        function Mu = applyMdown(obj, uBdy)
            % function that applies M, the integral matrix in matrix free form, it
            % only relies on local solves to apply the matrix.
            
             
            Mu = zeros(size(uBdy));
            
            ii = 1;
            % parsing the information at the boundary
            indlocal1 = (1:obj.M{ii}.n1);  % defining the local indeces
            indlocaln = (1:obj.M{ii}.nn);
            indStart = 0; % defining the last index
            indShift1 = obj.M{ii}.n1;
            indShiftn = obj.M{ii}.nn;
            % parsing the vector uBdy
            uy0  = 0*uBdy(indlocal1);
            uy1  = 0*uBdy(indlocal1);
            uyN  = uBdy(          indlocaln);
            uyNp = uBdy(indShiftn+indlocaln);

            %[vy0, vy1, vyN, vyNp] = obj.M{1}.applyBlockOperator( uy0, uy1, uyN, uyNp );
            Mu(indlocaln) = - uyN;
            
            for ii = 2:obj.nLayer-1
                % parsing the information at the boundary
                indlocal1 = indStart + (1:obj.M{ii}.n1);  % defining the local indeces
                indlocaln = indStart + 2*obj.M{ii}.n1 + (1:obj.M{ii}.nn);
                indStart =  max(indlocal1) + obj.M{ii}.n1;
                indShift1 = obj.M{ii}.n1;
                indShiftn = obj.M{ii}.nn;
                
               
                % parsing the vector uBdy
                uy0  = uBdy(          indlocal1);
                uy1  = uBdy(indShift1+indlocal1);
                uyN  = uBdy(          indlocaln);
                uyNp = uBdy(indShiftn+indlocaln);
                                
                [vy0, vy1, vyN, vyNp] = obj.M{ii}.applyBlockOperator( uy0, uy1, uyN, uyNp );
                Mu(indShift1+indlocal1) = vy1 - uy1;
                
                [vy0, vy1, vyN, vyNp] = obj.M{ii}.applyBlockOperator( uy0, uy1, 0*uyN, 0*uyNp );
                Mu(          indlocaln) = vyN - uyN;
                
            end
            
            ii = obj.nLayer;
            indlocal1 = indStart + (1:obj.M{ii}.n1);  % defining the local indeces
            indlocaln = indStart + 2*obj.M{ii}.n1 + (1:obj.M{ii}.nn);
            indStart =  max(indlocal1) + obj.M{ii}.n1;
            indShift1 = obj.M{ii}.n1;
            indShiftn = obj.M{ii}.nn;
            % parsing the vector
            uy0  = uBdy(         indlocal1);
            uy1  = uBdy(indShift1+indlocal1);
            uyN  = 0*uBdy(indlocal1);
            uyNp = 0*uBdy(indlocal1);
            
            [vy0, vy1, vyN, vyNp] = obj.M{ii}.applyBlockOperator( uy0, uy1, 0*uyN, 0*uyNp );
            Mu(indShift1+indlocal1)          = vy1 - uy1;
            
        end % function applyMdown
        
        function Mu = applyMup(obj, uBdy)
            % function that applies M, the integral matrix in matrix free form, it
            % only relies on local solves to apply the matrix.
            
            Mu = zeros(size(uBdy));
            
            ii = 1;
            % parsing the information at the boundary
            indlocal1 = (1:obj.M{ii}.n1);  % defining the local indeces
            indlocaln = (1:obj.M{ii}.nn);
            indStart = 0; % defining the last index
            indShift1 = obj.M{ii}.n1;
            indShiftn = obj.M{ii}.nn;
            % parsing the vector uBdy
            uy0  = 0*uBdy(indlocal1);
            uy1  = 0*uBdy(indlocal1);
            uyN  = uBdy(          indlocaln);
            uyNp = uBdy(indShiftn+indlocaln);

            [~, vy1, vyN, ~] = obj.M{1}.applyBlockOperator( uy0, uy1, uyN, uyNp );
            Mu(indlocaln) = vyN - uyN;
            
            for ii = 2:obj.nLayer-1
                % parsing the information at the boundary
                indlocal1 = indStart + (1:obj.M{ii}.n1);  % defining the local indeces
                indlocaln = indStart + 2*obj.M{ii}.n1 + (1:obj.M{ii}.nn);
                indStart =  max(indlocal1) + obj.M{ii}.n1;
                indShift1 = obj.M{ii}.n1;
                indShiftn = obj.M{ii}.nn;
                
               
                % parsing the vector uBdy
                uy0  = uBdy(          indlocal1);
                uy1  = uBdy(indShift1+indlocal1);
                uyN  = uBdy(          indlocaln);
                uyNp = uBdy(indShiftn+indlocaln);
                
                [~, vy1, vyN, ~]   = obj.M{ii}.applyBlockOperator( 0*uy0, 0*uy1, uyN, uyNp );
                Mu(indShift1+indlocal1) = vy1 - uy1;

                [~, vy1, vyN, ~]   = obj.M{ii}.applyBlockOperator( uy0, uy1, uyN, uyNp );
                Mu(          indlocaln) = vyN - uyN;
                
            end
            
            ii = obj.nLayer;
            indlocal1 = indStart + (1:obj.M{ii}.n1);  % defining the local indeces
            indlocaln = indStart + 2*obj.M{ii}.n1 + (1:obj.M{ii}.nn);
            indStart =  max(indlocal1) + obj.M{ii}.n1;
            indShift1 = obj.M{ii}.n1;
            indShiftn = obj.M{ii}.nn;
            % parsing the vector
            uy0  = uBdy(         indlocal1);
            uy1  = uBdy(indShift1+indlocal1);
            uyN  = 0*uBdy(indlocal1);
            uyNp = 0*uBdy(indlocal1);
            
            Mu(indShift1+indlocal1)   = - uy1;
       
        end % function applyMup
        
        function Lu = applyL(obj, uBdy)
            % function that applies L, the lower triangular par of the polired system
            % after a suitable permutation of the unknowns 
            
            Lu = zeros(size(uBdy));

            
           
            indStart = 0;
            
            
            for ii = 2:obj.nLayer-1
                % parsing the information at the boundary
                indlocal1 = indStart + (1:obj.M{ii}.n1);  % defining the local indeces
                indlocaln = indStart + 2*obj.M{ii}.n1 + (1:obj.M{ii}.nn);
                indStart =  max(indlocal1) + obj.M{ii}.n1;
                indShift1 = obj.M{ii}.n1;
                indShiftn = obj.M{ii}.nn;
                
               
                % parsing the vector uBdy
                uy0  = uBdy(          indlocal1);
                uy1  = uBdy(indShift1+indlocal1);
                uyN  = uBdy(          indlocaln);
                uyNp = uBdy(indShiftn+indlocaln);
                
                [vy0, vy1, vyN, vyNp]   = obj.M{ii}.applyBlockOperator( uy0, uy1, uyN, uyNp );
                Lu(         indlocal1) = vy0;
                Lu(indShift1+indlocal1) = vy1 - uy1;
                
            end
            
            ii = obj.nLayer;
            indlocal1 = indStart + (1:obj.M{ii}.n1);  % defining the local indeces
            indlocaln = indStart + 2*obj.M{ii}.n1 + (1:obj.M{ii}.nn);
            indStart =  max(indlocal1) + obj.M{ii}.n1;
            indShift1 = obj.M{ii}.n1;
            indShiftn = obj.M{ii}.nn;
            % parsing the vector
            uy0  = uBdy(         indlocal1);
            uy1  = uBdy(indShift1+indlocal1);
            uyN  = 0*uBdy(indlocal1);
            uyNp = 0*uBdy(indlocal1);
            
            [vy0, vy1, vyN, vyNp] = obj.M{ii}.applyBlockOperator( uy0, uy1, uyN, uyNp );
            Lu(         indlocal1)  = vy0;
            Lu(indShift1+indlocal1) = vy1 - uy1;

        end % function applyL


        function Mu = applyM0down(obj,uBdy)
            % function that applies M0 down, the integral matrix in matrix free form, it
            % only relies on local solves to apply the matrix.
               
            Mu = zeros(size(uBdy));
            
            ii = 1;
            % parsing the information at the boundary
            indlocal1 = (1:obj.M{ii}.n1);  % defining the local indeces
            indlocaln = (1:obj.M{ii}.nn);
            indStart = 0; % defining the last index
            indShift1 = obj.M{ii}.n1;
            indShiftn = obj.M{ii}.nn;
            % parsing the vector uBdy
            uy0  = 0*uBdy(indlocal1);
            uy1  = 0*uBdy(indlocal1);
            uyN  = uBdy(          indlocaln);
            uyNp = uBdy(indShiftn+indlocaln);

            %[vy0, vy1, vyN, vyNp] = obj.M{1}.applyBlockOperator( uy0, uy1, uyN, uyNp );
            Mu(indlocaln) = - uyNp;
            
            for ii = 2:obj.nLayer-1
                % parsing the information at the boundary
                indlocal1 = indStart + (1:obj.M{ii}.n1);  % defining the local indeces
                indlocaln = indStart + 2*obj.M{ii}.n1 + (1:obj.M{ii}.nn);
                indStart =  max(indlocal1) + obj.M{ii}.n1;
                indShift1 = obj.M{ii}.n1;
                indShiftn = obj.M{ii}.nn;
                
               
                % parsing the vector uBdy
                uy0  = uBdy(          indlocal1);
                uy1  = uBdy(indShift1+indlocal1);
                uyN  = uBdy(          indlocaln);
                uyNp = uBdy(indShiftn+indlocaln);
                
                [vy0, vy1, vyN, vyNp] = obj.M{ii}.applyBlockOperator( uy0, uy1, uyN, uyNp );
                Mu(indShift1+indlocal1) = vy0;
                
                [vy0, vy1, vyN, vyNp] = obj.M{ii}.applyBlockOperator( uy0, uy1, 0*uyN, 0*uyNp );
                Mu(          indlocaln) = vyNp - uyNp;
                
                
            end
            
            ii = obj.nLayer;
            indlocal1 = indStart + (1:obj.M{ii}.n1);  % defining the local indeces
            indlocaln = indStart + 2*obj.M{ii}.n1 + (1:obj.M{ii}.nn);
            indStart =  max(indlocal1) + obj.M{ii}.n1;
            indShift1 = obj.M{ii}.n1;
            indShiftn = obj.M{ii}.nn;
            % parsing the vector
            uy0  = uBdy(         indlocal1);
            uy1  = uBdy(indShift1+indlocal1);
            uyN  = 0*uBdy(indlocal1);
            uyNp = 0*uBdy(indlocal1);
            
            [vy0, vy1, vyN, vyNp] = obj.M{ii}.applyBlockOperator( uy0, uy1, uyN, uyNp );
            Mu(indShift1+indlocal1)   = vy0;
            
        end % function applyM0down 
          
        function Mu = applyM0up(obj, uBdy)
            % function that applies M, the integral matrix in matrix free form, it
            % only relies on local solves to apply the matrix.
            Mu = zeros(size(uBdy));
            
            ii = 1;
            % parsing the information at the boundary
            indlocal1 = (1:obj.M{ii}.n1);  % defining the local indeces
            indlocaln = (1:obj.M{ii}.nn);
            indStart = 0; % defining the last index
            indShift1 = obj.M{ii}.n1;
            indShiftn = obj.M{ii}.nn;
            % parsing the vector uBdy
            uy0  = 0*uBdy(indlocal1);
            uy1  = 0*uBdy(indlocal1);
            uyN  = uBdy(          indlocaln);
            uyNp = uBdy(indShiftn+indlocaln);

            [vy0, vy1, vyN, vyNp] = obj.M{1}.applyBlockOperator( uy0, uy1, uyN, uyNp );
            Mu(indlocaln) = vyNp;
            
            for ii = 2:obj.nLayer-1
                % parsing the information at the boundary
                indlocal1 = indStart + (1:obj.M{ii}.n1);  % defining the local indeces
                indlocaln = indStart + 2*obj.M{ii}.n1 + (1:obj.M{ii}.nn);
                indStart =  max(indlocal1) + obj.M{ii}.n1;
                indShift1 = obj.M{ii}.n1;
                indShiftn = obj.M{ii}.nn;
                
               
                % parsing the vector uBdy
                uy0  = uBdy(          indlocal1);
                uy1  = uBdy(indShift1+indlocal1);
                uyN  = uBdy(          indlocaln);
                uyNp = uBdy(indShiftn+indlocaln);
                
                [vy0, vy1, vyN, vyNp] = obj.M{ii}.applyBlockOperator( 0*uy0, 0*uy1, uyN, uyNp );
                Mu(indShift1+indlocal1) = vy0 - uy0;
                
                [vy0, vy1, vyN, vyNp] = obj.M{ii}.applyBlockOperator( uy0, uy1, uyN, uyNp );
                Mu(          indlocaln) = vyNp;
                
            end
            
            ii = obj.nLayer;
            indlocal1 = indStart + (1:obj.M{ii}.n1);  % defining the local indeces
            indlocaln = indStart + 2*obj.M{ii}.n1 + (1:obj.M{ii}.nn);
            indStart =  max(indlocal1) + obj.M{ii}.n1;
            indShift1 = obj.M{ii}.n1;
            indShiftn = obj.M{ii}.nn;
            % parsing the vector
            uy0  = uBdy(         indlocal1);
            uy1  = uBdy(indShift1+indlocal1);
            uyN  = 0*uBdy(indlocal1);
            uyNp = 0*uBdy(indlocal1);
            
            Mu(indShift1+indlocal1)  = -uy0;
  
        end % function applyM0up
        
        function Mu = applyDdownInv(obj, uBdy)
            % function that inverse of Ddown (without permutation)
            % to uBdy
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% WE NEED TO BE CAREFULL IN HERE!!
            
            Mu = zeros(size(uBdy));
            
            ii = 1;
            % parsing the information at the boundary
            indlocal1 = (1:obj.M{ii}.n1);  % defining the local indeces
            indlocaln = (1:obj.M{ii}.nn);
            indStart = 0; % defining the last index
            indShift1 = obj.M{ii}.n1;
            indShiftn = obj.M{ii}.nn;
            % parsing the vector uBdy
            uyN  = uBdy(          indlocaln);
            uyNp = uBdy(indShiftn+indlocaln);
            
            %[~, ~, vyN, vyNp] = obj.M{1}.applyBlockOperator( uy0, uy1, uyN, uyNp );
            Mu(          indlocaln) = - uyN;
            Mu(indShiftn+indlocaln) = - uyNp;
            
            for ii = 2:obj.nLayer-1
                % parsing the information at the boundary
                indlocal1 = indStart + (1:obj.M{ii}.n1);  % defining the local indeces
                indlocaln = indStart + 2*obj.M{ii}.n1 + (1:obj.M{ii}.nn);
                indStart =  max(indlocal1) + obj.M{ii}.n1;
                indShift1 = obj.M{ii}.n1;
                indShiftn = obj.M{ii}.nn;
                
                % parsing the vector uBdy
                vyN  = Mu(            indlocal1);
                vyNp = Mu(indShift1 + indlocal1);
                uyN  = uBdy(         indlocaln);
                uyNp = uBdy(indShiftn+indlocaln);
                [~, ~ , vyN, vyNp] = obj.M{ii}.applyBlockOperator( vyN, vyNp, 0*uyN, 0*uyNp );
                Mu(         indlocaln) = vyN  - uyN;
                Mu(indShiftn+indlocaln) = vyNp - uyNp;
                
            end

        end % function applyDdownInv
        
        function Mu = applyDupInv(obj, uBdy)
            % function that applies M, the integral matrix in matrix free form, it
            % only relies on local solves to apply the matrix.
            
            Mu = zeros(size(uBdy));
            
            ii = obj.nLayer;
            % parsing the information at the boundary
            indStart = length(uBdy)-2*obj.M{ii}.n1; % defining the first index
            indlocal = indStart + (1:obj.M{ii}.n1);  % defining the local indeces
            indShift = obj.M{ii}.n1;
            
            % parsing the vector uBdy
            uy0  = uBdy(indlocal);
            uy1  = uBdy(indShift+indlocal);

            %[vy0, vy1, vyN, vyNp] = obj.M{1}.applyBlockOperator( uy0, uy1, uyN, uyNp );
            Mu(indlocal)             = - uy0;
            Mu(indShift+indlocal)    = - uy1;
            
            vy0 = -uy0;
            vy1 = -uy1;
                        
            for ii = obj.nLayer-1:-1:2
                % parsing the information at the boundary
                indlocal  = (indStart - 2*obj.M{ii}.nn )...
                             + (1:obj.M{ii}.nn);  % defining the local indeces
                indStart  = min(indlocal)-1; % carefull with the minus one here!!!
                indShift  = obj.M{ii}.nn;
                
                % parsing the vector uBdy
                uy0  = uBdy(           indlocal);
                uy1  = uBdy(  indShift+indlocal);
                
                [vy0, vy1, ~, ~ ] = obj.M{ii}.applyBlockOperator( 0*uy0, 0*uy1, vy0, vy1 );
                Mu(         indlocal) = vy0 - uy0; 
                Mu(indShift+indlocal) = vy1 - uy1;
                
                vy0 = Mu(         indlocal);
                vy1 = Mu(indShift+indlocal);
            end

        end % function applyDupInv
        
        function Dinvu = applyDinv(obj, uBdy, flag)
            % apply the inverse of of the diagonal of the polarized method

            uBdy = obj.P*uBdy;
            Dinvup   = obj.applyDupInv(uBdy(end/2+1:end));
            Dinvdown = obj.applyDdownInv(uBdy(1:end/2)); 
            Dinvu =  [Dinvdown; Dinvup];
            
        end %function applyDinv

        function Dinvu = applyGSinv(obj, uBdy, flag)
            % perform one iteration of the Gauss-Seidel iteration 
            % as a preconditioner

            uBdy = obj.P*uBdy;
            Dinvdown = obj.applyDdownInv(uBdy(1:end/2)); 
            aux      = obj.applyL(Dinvdown);
            Dinvup   = obj.applyDupInv(uBdy(end/2+1:end)- aux);
            Dinvu =  [Dinvdown; Dinvup];

        end %function applyDinv
        
        
        
        function MMu = applyMpolarized(obj, uBdy)
            fprintf('Applying the polarized matrix within the Gmres iteration \n')
           % apply the polarized system 
           uBdyDown = uBdy(1:end/2);
           uBdyUp   = uBdy(end/2+1:end);
           
           MMu = [obj.applyMdown(uBdyDown)  + obj.applyMup(uBdyUp) ;
                  obj.applyM0down(uBdyDown) + obj.applyM0up(uBdyUp)];
            
            
        end
        
        
        
        function compressGreenFunctions(obj, max_rank, epsilon)
            % function to compress all the Green's Functions
            fprintf('Compressing the Green functions \n')
           for ii =1:obj.nLayer
               obj.M{ii}.compressGreenFunctions(max_rank, epsilon);
           end
        end
        
        % This needs to be modified in order to consider the different
        % sizes 
        function initPermutationMatrix(obj)
            % initializes the permutation matrix
            nLayer = obj.nLayer;
            E = speye(4*(nLayer-1));
            p_aux   = kron(1:2:2*(nLayer-1)-1, [1 1]) + kron(ones(1,nLayer-1), [0 2*(nLayer-1) ]);
            p_aux_2 = kron(2:2:2*(nLayer-1), [1 1])   + kron(ones(1,nLayer-1), [2*(nLayer-1) 0]);
            p = E( [p_aux, p_aux_2],: );
            obj.P = kron(p, speye(obj.M{1}.nn));
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
                for ii = 1:1: obj.nLayer
                    s_array{ii} = obj.M{ii}.extractLocalSource(s);
                end
            elseif nargin == 3
                if strcmp(flag,'extended') % for the single partition test
                    s_array = {};
                    for ii = 1:1: obj.nLayer
                        % we have to split the source and add one more pixel up
                        % and down
                        s_aux = s(find( (obj.Y>= obj.y_sample(ii)-obj.h-10*eps).*(obj.Y<=obj.x_extrp(ii)+obj.h+10*eps)));
                        s_aux = reshape(s_aux, numel(s_aux)/obj.nn, obj.nn);
                        % padding zeros up and down (this can be a problem if the
                        % source is no completely inside the domain
                        s_array{ii} = [ zeros(obj.npml-1,obj.nn); s_aux ; zeros(obj.npml-1,obj.nn) ];
                    end
                end
            end
                
            end % function source_splitting
            
        
        
        % Function to concatenate a series of partial local solution in
        % order to obtain a global solution. Its main application is to
        % compare different solutions 
        
        
        function uGlobal = reconstruct(obj, uArray)
            % fucntion to reconstruct the global wavefield from the local
            % reconstructions 
            uGlobal = zeros(size(obj.MBase.node(obj.MBase.freeNode,:),1),1);
            for ii = 1:obj.nLayer;
               uGlobal(obj.M{ii}.indIntGlobal) = uArray{ii}(obj.M{ii}.indIntLocal);
            end
            
        end
        
        
%         function u_global = concatenate(obj, u_array, flag)
%             
%             aux_counter = 0;
%             if nargin==2
%                 u_global = zeros(size(obj.Xi));
%                 for ii = 1:numel(u_array)
%                     % we should have continuity in the overlaping interfacing
%                     % so this should be more than enough. Otherwise I need to
%                     % change this in order to consider the jump conditions.
%                     
%                     u_reordered = reshape(u_array{ii}, obj.M{ii}.ny, obj.M{ii}.nn);
%                     u_global( aux_counter+1:obj.M{ii}.nyi+aux_counter, 1:obj.M{ii}.nni) = ...
%                         u_reordered( obj.M{ii}.npml+1:obj.M{ii}.npml+obj.M{ii}.nyi, obj.M{ii}.npml+1:obj.M{ii}.npml+obj.M{ii}.nni );
%                     aux_counter = aux_counter + (obj.M{ii}.nyi) ;
%                 end
%                 
%             elseif nargin==3 && strcmp(flag,'extended') % this flag is to add the PML 
%                 u_global = zeros(size(obj.X));
%                 for ii = 1:numel(u_array)
%                     % we should have continuity in the loverlaping interfacing
%                     % so this should be more than enough. Otherwise I need to
%                     % change this in order to consider the jump conditions.
%                     
%                     % have to reorder the solutions (this is legacy from
%                     % Vincent;s code)
%                     u_reordered = reshape(u_array{ii}, obj.M{ii}.ny, obj.M{ii}.nn);
%                     if ii~= 1 && ii~= numel(u_array)
%                         
%                         u_global( aux_counter+1:obj.M{ii}.nyi+aux_counter, 1:obj.M{ii}.nn) = ...
%                             u_reordered( obj.M{ii}.npml+1:obj.M{ii}.npml+obj.M{ii}.nyi, 1:obj.M{ii}.nn);
%                         aux_counter = aux_counter + obj.M{ii}.nyi  ;
%                         
%                     elseif ii == 1 % in the first one we need to add the pml at the begining
%                         
%                         u_global( aux_counter+1:obj.M{ii}.nyi+aux_counter+obj.M{ii}.npml, 1:obj.M{ii}.nn) = ...
%                             u_reordered(1:obj.M{ii}.npml+obj.M{ii}.nyi, 1:obj.M{ii}.nn);
%                         aux_counter = aux_counter + obj.M{ii}.npml + obj.M{ii}.nyi ;
%                         
%                     elseif ii == numel(u_array) % last one, have to add the pml as well
%                         
%                         u_global( aux_counter+1:obj.M{ii}.nyi+aux_counter+obj.M{ii}.npml, 1:obj.M{ii}.nn) = ...
%                             u_reordered( obj.M{ii}.npml+1:2*obj.M{ii}.npml+obj.M{ii}.nyi, 1:obj.M{ii}.nn);
%                         aux_counter = aux_counter + obj.M{ii}.nyi ;
%                         
%                     end
%                 end
%                 
%                 
%             end
%             
%             
%             
%         end % function concatenate
%         
        %^^^^^^^^^^^^^^^^^^^^^^^ TO DO
        
         function uArraySol = solveLocal(obj,s,uBdy)
            % function reconstruct the exact solution locally once the global
            % traces have already been computed.
            
            ii = 1;

            indlocaln = (1:obj.M{ii}.nn);
            indStart = 0; % defining the last index
            indShiftn = obj.M{ii}.nn;
            
            s_array = obj.source_splitting( s);
            
            uy0  = zeros(indShiftn,1);
            uy1  = uy0; 
            indlocaln = (1:obj.M{ii}.nn);
            uyN  = uBdy(indlocaln);
            uyNp = uBdy( indShiftn + indlocaln);
            
            u_aux       = obj.M{ii}.solveTraces(s_array{ii}, uy0, uy1, uyN, uyNp);
            uArraySol{ii} = u_aux;
            
            for ii = 2:obj.nLayer-1
               
                indlocal1 = indStart + (1:obj.M{ii}.n1);  % defining the local indeces
                indlocaln = indStart + 2*obj.M{ii}.n1 + (1:obj.M{ii}.nn);
                indStart =  max(indlocal1) + obj.M{ii}.n1;
                indShift1 = obj.M{ii}.n1;
                indShiftn = obj.M{ii}.nn;
                
               
                % parsing the vector uBdy
                uy0  = uBdy(          indlocal1);
                uy1  = uBdy(indShift1+indlocal1);
                uyN  = uBdy(          indlocaln);
                uyNp = uBdy(indShiftn+indlocaln);
                
                % solve locallyusing the forcing
                u_aux       = obj.M{ii}.solveTraces(s_array{ii}, uy0, uy1, uyN, uyNp);
                uArraySol{ii} = u_aux;
            end
            
            ii = obj.nLayer;
            
            indlocal1 = indStart + (1:obj.M{ii}.n1);  % defining the local indeces
            indShift1 = obj.M{ii}.n1;
            
            % parsing the vector
            uy0  = uBdy(          indlocal1);
            uy1  = uBdy(indShift1+indlocal1);
            uyN  = 0*uBdy(indlocal1);
            uyNp = 0*uBdy(indlocal1);
            
            
            u_aux       = obj.M{ii}.solveTraces(s_array{ii}, uy0, uy1, uyN, uyNp);
            uArraySol{ii} = u_aux;
            
            
        end
        
        %^^^^^^^^^^^^^^^^^^^^^^^ TO DO
               
       %TODO function to compute the Green's function inside each
       %horizontal layer
        
        % function to display the local wavefields given by solving the
        % local Systems
        function DisplayField(obj, u_array)
            % have to check the lenght of an array 
            for ii = 1:obj.nLayer
                figure(100+ii);
                DisplayField(u_array{ii}, obj.M{ii}.x, obj.M{ii}.y, obj.M{ii}.npml);
            end
        end % function DisplayField
        
        %small function that converts the sparse data on a 
        function s = struct2field(obj, s_struct)
            s = zeros(size(obj.X));
            for ii = 1:numel(s_struct.y)
               ind_y = find(abs(obj.Y-s_struct.y(ii))<2*eps);
               s(ind_y) = s_struct.u(:,ii);
            end
            
        end
            
    end % methods
    
end % classdef

 