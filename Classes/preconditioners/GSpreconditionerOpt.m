function uPrecond = GSpreconditionerOpt(MArray, indx0, indx1, indxn ,indxnp, ...
                                     indIntGlobal, indIntLocal, fInt )
% function uPrecond = GSpreconditioner(MArray, indx0, indx1, indxn ,indxnp,
%                                      indIntGlobal, indIntLocal, fInt )
% This function applies the Polarized traces preconditioner using the
% optimized Gauss-Seidel version. (See Fast Domain Decomposition solver for
% the Lippmann-Schwinger equation). 
% TODO: this needs to be futher encapsulated                               
                                 
%% extracting some usefull information   
nSub =  length(MArray);

% size of the local traces
n = size(indxn{1},1);

u_0  = zeros(n*nSub,1);
u_1  = zeros(n*nSub,1);
u_n  = zeros(n*nSub,1);
u_np = zeros(n*nSub,1);

index = 1:n;

%% source partitionning 

fIntLocal = {};

% we need to extract the sizes of each local solution 
localSizes = zeros( 1, nSub);


%% local solve + downwards sweep
for ii = 1:nSub
   
    % making a local copy of the local rhs
    fIntLocal{ii} = zeros(size(MArray{ii}.node(MArray{ii}.freeNode,1),1),1);
    fIntLocal{ii}(indIntLocal{ii}) = fInt(indIntGlobal{ii});
    localSizes(ii) = size(indIntGlobal{ii},1);

    if ii ~=1
        fIntLocal{ii}(indx0{ii}) =  fIntLocal{ii}(indx0{ii}) +...
                   MArray{ii}.H(indx0{ii},indx1{ii})*u_np((ii-2)*n + index);
        fIntLocal{ii}(indx1{ii}) = fIntLocal{ii}(indx1{ii}) -...
                    MArray{ii}.H(indx1{ii},indx0{ii})*u_n((ii-2)*n + index);
    end

        % solving the rhs
        vDown = MArray{ii}.solveInt(fIntLocal{ii});

        % extracting the traces
        if ii ~= nSub
            u_n((ii-1)*n  + index) = vDown(indxn{ii});
            u_np((ii-1)*n + index) = vDown(indxnp{ii});
        end
end

localLim = [0 cumsum(localSizes)];

%% upwards sweep + reflections + reconstruction

uPrecond = zeros(size(fInt,1),1);

     for ii = nSub:-1:1
      
        if ii~= nSub
            fIntLocal{ii}(indxnp{ii}) = fIntLocal{ii}(indxnp{ii}) + ...
                MArray{ii}.H(indxnp{ii},indxn{ii})*u_0((ii)*n + index);
            fIntLocal{ii}(indxn{ii})  = fIntLocal{ii}(indxn{ii}) - ...
                MArray{ii}.H(indxn{ii},indxnp{ii})*u_1((ii)*n + index);
        end
 
        % solving the local problem
        uUp = MArray{ii}.solveInt(fIntLocal{ii});

        if ii > 1
            u_0((ii-1)*n + index) = uUp(indx0{ii});
            u_1((ii-1)*n + index) = uUp(indx1{ii}) - u_np((ii-2)*n + index);
        end

        % reconstructing the problem on the fly
        uPrecond(localLim(ii)+1:localLim(ii+1)) = uUp(indIntLocal{ii});
     end
    
end