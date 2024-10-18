function P = computeCG1toDG1(mesh)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$
    % computes prolongation from CG1 to DG1
    %
    % Input:
    %     mesh:  current mesh 
    %   
    % Output: 
    %        P:  prolongation matrix
    %
    % M. Hauck, A. Lozinski
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$

    [np,d] = size(mesh.p);
    nt = size(mesh.t,1);
    P = sparse(1:nt*(d+1),reshape(mesh.t.',1,[]),1,nt*(d+1),np);
end % function