function mesh = getMesh(geom)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get initial mesh
    %
    % Input:
    %     geom:  geometry of mesh to load 
    %   
    % Output: 
    %     mesh:  initial mesh
    %
    % M. Hauck, A. Lozinski
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$
    
    switch geom
        case 'Interval'
            p = [0;1];
            t = [1 2];
        case 'Triangle'
            p = [0 0; 1 0; 1 1];
            t = [1 2 3];
        case 'Square'
            p = [0 0; 1 0; 1 1; 0 1];
            t = [1 2 3; 1 3 4];
        case 'Cube'
            p = [0 0 0; 0 1 0; 1 0 0; 1 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 1];
            t = [6 2 3 1; 6 3 5 1; 6 2 4 3; 6 4 8 3; 6 7 5 3; 6 8 7 3];
        otherwise
            error('geometry not found.');
    end % switch

    % edge data 
    [e,te,e2t,p2e] = getEdges(struct('p',p,'t',t));

    % create mesh
    mesh = struct('p',p, ...            % coordinates
                  't',t,...             & node indices of elements
                  'e',e,...             % node indices of edges
                  'np',size(p,1),...    % number of nodes
                  'nt',size(t,1),...    % number of elements 
                  'ne',size(e,1),...    $ number of edges
                  'te',te,...           % edge indices of elements
                  'e2t',e2t,...         % in which elements is an edge
                  'p2e',p2e);           % in which edge is a node
end % function