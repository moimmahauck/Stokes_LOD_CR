function c = getPatches(mesh,ell)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get incidence matrix for ell-th order element patches
    %
    % Input:
    %     mesh:  mesh of which we want to compute the patches
    %      ell:  oversampling parameter
    %   
    % Output: 
    %        c:  incidence matrix C with C(i,j) = 1 iff intersection of 
    %            elements T_j and T_k is not empty
    %
    % M. Hauck, A. Lozinski
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$

    % compute element patches
    nt = size(mesh.t,1);
    [np,d] = size(mesh.p);
    C = sparse(mesh.t,repmat((1:nt)',1,d+1),1,np,nt);
    C = spones(C'*C); 
    % patches are encoded by matrix c \in {true,false}^NTHxNTH: 
    % c(j,k) = 1 <=> T_j subset N^nloc(T_k) and T_j not equal to T_k
    c = speye(nt);
    for k = 1:ell
        c = C*c;
    end % for
    % remove seed element, added back later
    c = logical(spones(c)-speye(size(c))); 
end % function