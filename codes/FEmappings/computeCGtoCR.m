function P = computeCGtoCR(mesh)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$
    % computes prolongation from CG1 to CR
    %
    % Input:
    %     mesh:  current mesh 
    %   
    % Output: 
    %        P:  prolongation matrix
    %
    % M. Hauck, A. Lozinski
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$

    P = sparse(repmat(1:mesh.ne,2,1)',mesh.e,.5*ones(size(mesh.e)),mesh.ne,mesh.np);
end % function