function mids = computeElemMids(mesh)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute midpoint of elements of given mesh
    %
    % Input:
    %     mesh:  mesh 
    %   
    % Output: 
    %     mids:  midpoints of elements
    %
    % M. Hauck, A. Lozinski
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$

    mids = (mesh.p(mesh.t(:,1),:) + mesh.p(mesh.t(:,2),:) + mesh.p(mesh.t(:,3),:))/3;
end % function