function S = assembleStiffness_CR(mesh,coeff)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$
    % assembles possibly weighted CR stiffness matrix
    %
    % Input:
    %     mesh:  simplicial mesh with points p and elements t
    %    coeff:  piecewise constant scalar coefficient vector
    %
    % Output: 
    %        S:  stiffness matrix
    %
    % M. Hauck, A. Lozinski
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$

    % set default coefficient to 1
    if nargin < 2 || isempty(coeff)
        coeff = ones(size(mesh.t,1),1);
    end % if
    
    G = computeGrads_CR(mesh);

    % assemble local stiffness matrices
    Sloc = zeros(mesh.nt,9);
    c = coeff.*simpvol(mesh);
    
    for k1 = 1:3
        for k2 = 1:3
            Sloc(:,k2+3*(k1-1)) = c.*dot(G(:,:,k1),G(:,:,k2),2);
        end % for
    end % for

    % assemble global matrix
    I = (1:3)'*ones(1,3);
    S = sparse(mesh.te(:,reshape(I',1,9)),...
               mesh.te(:,reshape(I,1,9)),Sloc,mesh.ne,mesh.ne);
end % function