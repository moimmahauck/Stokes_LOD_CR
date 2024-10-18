function P = computeDG0toDG1(mesh)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$
    % computes prolongation from DG0 to DG1
    %
    % Input:
    %     mesh:  current mesh 
    %   
    % Output: 
    %        P:  prolongation matrix
    %
    % M. Hauck, A. Lozinski
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$

    d = size(mesh.p,2);
    P = kron(speye(mesh.nt),ones(d+1,1));   
end % function