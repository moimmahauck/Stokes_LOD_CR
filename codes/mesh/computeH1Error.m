function errH1 = computeH1Error(mesh,coeff,Dfun)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute H1-error between discrete and continuous function
    %
    % Input:
    %     mesh:  mesh 
    %     coeff: coefficients of discrete function
    %     Dfun:  derivative of continuous function as function handle 
    %   
    % Output: 
    %    errH1:  value of H1-error
    %
    % M. Hauck, A. Lozinski
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$

    coeff = reshape(coeff,3,[])';
    
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
    
    % Assemble local mass matrices: (u,v)_T->int_T rho1*rho2*u*v
    errH1loc = zeros(mesh.nt,2);
    vol = simpvol(mesh);

    % loop over local pairs of dofs (parallel in t)
    valatquad = computeGrads_P1(mesh,coeff);
    for q = 1:length(W)
        xQuad = (mesh.p(mesh.t(:,2),:)-mesh.p(mesh.t(:,1),:)).*X(q,1) + ...
                (mesh.p(mesh.t(:,3),:)-mesh.p(mesh.t(:,1),:)).*X(q,2) + ...
                 mesh.p(mesh.t(:,1),:);
        
        tmp = Dfun(xQuad);
        for j = 1:2
            errH1loc(:,j) = errH1loc(:,j) + W(q)*vol.*(valatquad(:,j)-tmp(:,j)).^2;
        end % for
    end % for
    errH1 = sqrt(sum(errH1loc(:)));
end % function
