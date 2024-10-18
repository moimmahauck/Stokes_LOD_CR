function [meshh,P,P0,P1dg] = refineMesh(meshH,nref)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % refines a coarse triangulation TH nref times by applying the function
    % refine nref times. Suitable restriction and prolongation operator
    % between corresponding mesh functions are computed
    %
    % Input:
    %       TH:  coarse mesh (structure array consisting of array of vertices p
    %            and array of elements represented by vertex indices
    %     nref:  number of refinments [default = 1]
    % 
    % Output: 
    %       Th:  fine mesh after refinement
    %        P:  Prolongation matrix P1
    %            extrapolates P1 grid functions from the coarse to the fine mesh
    %       P0:  Prolongation matrix P0
    %            extrapolates P0 mesh functions from the coarse to the fine mesh
    %     P1dg:  Prolongation matrix P1dg
    %            extrapolates P1dg mesh functions from the coarse to the fine mesh
    %
    % M. Hauck, A. Lozinski
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$

    % Construct data structure
    if nargin < 2
        nref = 1;
    end % if
    
    % compute fine reference mesh
    P = speye(meshH.np);
    P0 = speye(meshH.nt);
    if nargout > 3
        P1dg = speye(meshH.nt*size(meshH.t,2));
    end % if
    
    meshh = meshH;

    for k = 1:nref
        if nargout > 3
            [meshh,p,p0,p1dg] = refine(meshh);
            P1dg = p1dg*P1dg;
        else
            [meshh,p,p0] = refine(meshh);
        end % if
        P = p*P;
        P0 = p0*P0;
    end % for
end % function