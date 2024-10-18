function errL2 = computeL2Error(mesh,coeff,fun)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute L2-error between discrete and continuous function
    %
    % Input:
    %     mesh:  mesh 
    %     coeff: coefficients of discrete function
    %      fun:  continuous function as function handle 
    %   
    % Output: 
    %    errL2:  value of L2-error
    %
    % M. Hauck, A. Lozinski
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$

    coeff = reshape(coeff,3,[])';
    
    % P1-FEM basis functions
    lam1 = @(X) 1-X(:,1)-X(:,2);
    lam2 = @(X) X(:,1);
    lam3 = @(X) X(:,2);

    evalbasisP1 = @(X) [lam1(X) lam2(X) lam3(X)];
    
    % quadrature rule
    c1 = [0 0]; c2 = [1 0]; c3 = [0 1];
    a1 = 0.445948490915965;
    a2 = 0.091576213509771;

    X = [a1*c1(1)+a1*c2(1)+(1-2*a1)*c3(1) a1*c1(2)+a1*c2(2)+(1-2*a1)*c3(2); ...
         a1*c1(1)+(1-2*a1)*c2(1)+a1*c3(1) a1*c1(2)+(1-2*a1)*c2(2)+a1*c3(2); ...
         (1-2*a1)*c1(1)+a1*c2(1)+a1*c3(1) (1-2*a1)*c1(2)+a1*c2(2)+a1*c3(2); ...
         a2*c1(1)+a2*c2(1)+(1-2*a2)*c3(1) a2*c1(2)+a2*c2(2)+(1-2*a2)*c3(2); ...
         a2*c1(1)+(1-2*a2)*c2(1)+a2*c3(1) a2*c1(2)+(1-2*a2)*c2(2)+a2*c3(2); ...
         (1-2*a2)*c1(1)+a2*c2(1)+a2*c3(1) (1-2*a2)*c1(2)+a2*c2(2)+a2*c3(2);];

    w1 = 0.223381589678010;
    w2 = 0.109951743655322;
    W = [w1; w1; w1; w2; w2; w2];

    % evaluate basis functions
    F = evalbasisP1(X).';

    % compute element integrals
    errL2loc = zeros(mesh.nt,1);
    vol = simpvol(mesh);

    for q = 1:length(W)
        xQuad = (mesh.p(mesh.t(:,2),:)-mesh.p(mesh.t(:,1),:)).*X(q,1) + ...
                (mesh.p(mesh.t(:,3),:)-mesh.p(mesh.t(:,1),:)).*X(q,2) + ...
                 mesh.p(mesh.t(:,1),:);
        
        valatquad = zeros(mesh.nt,1);
        for k = 1:3
            valatquad = valatquad + F(k,q)*coeff(:,k);
        end % for
        
        errL2loc = errL2loc + W(q)*vol.*(valatquad-fun(xQuad)).^2;
    end % for
    errL2 = sqrt(sum(errL2loc));
end % function
