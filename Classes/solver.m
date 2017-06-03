classdef solver < handle
    % wrapper for the LU factorization solver (based on UMFPACK)
    % with row and column pivoting and nested dissection reordering
    
    properties
        L   % L factor
        U   % U factor
        P   % row permutation matrix
        Q   % column permutation matrix
        n    % size of the matrix
        m
    end
    
    methods
        % so far we only consider the Lu factorization
        % we should makethis class an abstract class in order 
        % use the HF-ID 
        function obj = solver(H)
            % obj = solver(H)
            % function to instanciate a solver
            %
            % Input : H a matrix (with the boundary conditions
            %         already included inside
            %
            % Output: obj an instance of the solver class
            %         the object stores the LU factors and the permutation
            %         matrices such that P*H*Q = L*U
            
            fprintf('Factorizing the matrix using UMFPACK \n')
            % we perform the factorization and we save the LU factors
            [obj.L, obj.U, obj.P, obj.Q]=lu(H );    
            % we save the size of the matrix for debugging purpouses 
            [obj.n, obj.m] = size(H);
        end
        function u = solve(obj, f)
            % u = solve(obj, f)
            % function to solve Hu = f, using the direct sparse 
            % LU factorization algorihtm (UMFPACK)
            %
            % Input:  f  a matrix nxr containing r right-had sides
            %
            % Output: u the solution to the set of systems Hu = f

            
            % fprintf('Solving using LU factorization \n')
            if size(f,1) == obj.n  % if the dimension match
                u = obj.Q*(obj.U\(obj.L\(obj.P*f)));
            else
                fprintf('incorrect size of the right hand side \n');
            end
        end
    end
    
end
    
    