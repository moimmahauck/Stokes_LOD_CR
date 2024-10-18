function mids = computeEdgeMids(mesh)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute midpoint of edges of given mesh
    %
    % Input:
    %     mesh:  mesh 
    %   
    % Output: 
    %     mids:  midpoints of edges
    %
    % M. Hauck, A. Lozinski
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$
    
    mids = (mesh.p(mesh.e(:,1),:) + mesh.p(mesh.e(:,2),:))/2;
end % function