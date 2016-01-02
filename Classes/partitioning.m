classdef partitioning < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    % what are the properties that we need? 
    properties
    MArray
    indx0
    indx1
    indxn 
    indxnp
    indIntGlobal
    indIntLocal
    end
    
    methods
    function obj = partitioning(MArray, indx0, indx1, indxn ,indxnp, indIntGlobal, indIntLocal) 
            % obj = solver(H)
            % function to instanciate a solver
            %
            % Input : H a matrix (with the boundary conditions
            %         already included inside
            %
            % Output: obj an instance of the solver class
            %         the object stores the LU factors and the permutation
            %         matrices such that P*H*Q = L*U
            obj.MArray = MArray;
            obj.indx0 = indx0;
            obj.indx1 = indx1;
            obj.indxn = indxn; 
            obj.indxnp = indxnp;
            obj.indIntGlobal = indIntGlobal;
            obj.indIntLocal  = indIntLocal;
            
    end
    end
    
end