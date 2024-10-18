function M = assembleMass_CR(mesh,coeff)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$
    % assembles CR mass matrix
    %
    % Input:
    %     mesh:  simplicial mesh with points p and elements t
    %    coeff:  piecewise constant scalar coefficient vector
    %
    % Output: 
    %        M:  mass matrix
    %
    % M. Hauck, A. Lozinski
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$

    % set default coefficient to 1
    if nargin < 2 || isempty(coeff)
        coeff = ones(size(mesh.t,1),1);
    end % if
    
    % CR basis functions, quadrature points, and weights
    b = {@(x) 2*x(:,1)+2*x(:,2)-1, @(x) -2*x(:,1)+1, @(x) -2*x(:,2)+1};
    % this quadrature is exact for P2 functions -> no quadrature error
    X = [.5 0;.5 .5;0 .5];
    W = [1/3; 1/3; 1/3];
    
    % evaluate CR basis functions at quadrature points
    F = zeros(length(W),3);
    for j = 1:3
        F(:,j) = b{j}(X);
    end % for
    
    W = simpvol(mesh).*W';
    
    % assemble local mass matrices
    Mloc = zeros(mesh.nt,9);
    for k1 = 1:3
        for k2 = 1:3
            Mloc(:,k1+3*(k2-1)) = sum(W*(F(:,k1).*F(:,k2)),2);
        end % for
    end % for
    
    Mloc = coeff.*Mloc;

    % assemble global matrix
    I = (1:3)'*ones(1,3);
    M = sparse(mesh.te(:,reshape(I',1,9)),...
               mesh.te(:,reshape(I,1,9)),Mloc,mesh.ne,mesh.ne);
end % function