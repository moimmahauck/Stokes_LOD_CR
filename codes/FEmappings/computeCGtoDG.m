function P = computeCGtoDG(mesh,ord)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$
    % computes prolongation/projection from CG1 to DG0/DG1
    %
    % Input:
    %     mesh:  current mesh 
    %      ord:  order of DG space of projection/prolongation
    %   
    % Output: 
    %        P:  prolongation/projection matrix
    %
    % M. Hauck, A. Lozinski
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$

    switch ord
        case 0
            d = size(mesh.p,2);
            P = computeDG0toDG1(mesh).'*computeCG1toDG1(mesh)./(d+1);
        case 1
            P = computeCG1toDG1(mesh);
    end % switch
end % function