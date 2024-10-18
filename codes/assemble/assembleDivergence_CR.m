function D = assembleDivergence_CR(mesh)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$
    % assembles CR divergence matrix without quadrature error
    %
    % Input:
    %     mesh:  simplicial mesh with points p and elements t
    %   
    % Output: 
    %        D:  divergence matrix
    %
    % M. Hauck, A. Lozinski
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$

    G = computeGrads_CR(mesh);
    
    % assemble local matrices
    Dloc1 = zeros(mesh.nt,3);
    Dloc2 = zeros(mesh.nt,3);
    c = simpvol(mesh);
    
    for k = 1:3
        Dloc1(:,k) = c.*G(:,1,k);
        Dloc2(:,k) = c.*G(:,2,k);
    end % for
    
    % assemble global matrix
    D1 = sparse(repmat(1:mesh.nt,3,1)',mesh.te,Dloc1,mesh.nt,mesh.ne);
    D2 = sparse(repmat(1:mesh.nt,3,1)',mesh.te,Dloc2,mesh.nt,mesh.ne);
    D = [D1 D2];
end % function