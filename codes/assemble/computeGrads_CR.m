function G = computeGrads_CR(mesh,coeff)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$
    % computes piecewise gradient of CR grid function u 
    %
    % Input:
    %     mesh:  simplicial mesh
    %        u:  grid function
    %   
    % Output: 
    %        G: gradients (row-wise) according to the simplices in t
    %
    % M. Hauck, A. Lozinski
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$
    
    % midpoints of edges
    mids = computeEdgeMids(mesh);

    % normalized vectors connecting opposite nodes
    v1 = mids(mesh.te(:,2),:) - mids(mesh.te(:,3),:); v1 = v1./sqrt(sum(v1.^2,2));
    v2 = mids(mesh.te(:,3),:) - mids(mesh.te(:,1),:); v2 = v2./sqrt(sum(v2.^2,2));
    v3 = mids(mesh.te(:,1),:) - mids(mesh.te(:,2),:); v3 = v3./sqrt(sum(v3.^2,2));

    % height vectors of triangles from nodes
    h1 = mids(mesh.te(:,1),:) - mids(mesh.te(:,2),:); h1 = h1 - sum(v1.*h1,2).*v1;
    h2 = mids(mesh.te(:,2),:) - mids(mesh.te(:,1),:); h2 = h2 - sum(v2.*h2,2).*v2;
    h3 = mids(mesh.te(:,3),:) - mids(mesh.te(:,2),:); h3 = h3 - sum(v3.*h3,2).*v3;

    % gradients
    G = zeros(mesh.nt,2,3);
    G(:,:,1) = h1./sqrt(sum(h1.^2,2)).^2;
    G(:,:,2) = h2./sqrt(sum(h2.^2,2)).^2;
    G(:,:,3) = h3./sqrt(sum(h3.^2,2)).^2;

    if nargin == 2
        for k = 1:3
            G(:,:,k)=(coeff(:,k)*ones(1,2)).*G(:,:,k);
        end % for
        G = sum(G,3);
    end % if
end % function