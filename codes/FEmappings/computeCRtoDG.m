function P = computeCRtoDG(mesh)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$
    % computes prolongation from CR to DG
    %
    % Input:
    %     mesh:  current mesh 
    %   
    % Output: 
    %        P:  prolongation matrix
    %
    % M. Hauck, A. Lozinski
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$

    % maps to CR-DG representation
    P1 = sparse(reshape(1:numel(mesh.te),3,[])',mesh.te,ones(size(mesh.te)),3*mesh.nt,mesh.ne);
    % maps CG-DG representation to P1-DG
    loc = [-1 1 1;1 -1 1; 1 1 -1];
    P2 = kron(speye(mesh.nt),loc);
    % concatenating both operators
    P = P2*P1;
end % function