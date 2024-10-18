%% This code reproduces the first numerical experiment in the paper
%
% A Localized Orthogonal Decomposition Method for Heterogeneous Stokes
% Problems
%
% Note: the below parameter configuration was used on a cluster. The
% computation may not be feasuble on a laptop.
%
% M. Hauck, A. Lozinski

init;
close all; clear all;

ells = 1:6;
lvlHs = 3:5;

lvlh = 8;

%% create heterogeneous coefficient
alpha = .1;
beta  = 1;

Tc = refineMesh(getMesh('Square'),lvlh-2);
[Th,~,P0] = refineMesh(Tc,2);

mids = computeElemMids(Tc);

rng(27);
coeff = alpha + (beta-alpha)*rand(Tc.nt,1);
tf = abs(mids(:,2) - 3*mids(:,1).*(1-mids(:,1)))<2*2^-(lvlh-2) & mids(:,1) >.2;
coeff(tf) = 10;
coeff = P0*coeff;locerrH1 = zeros(length(lvlHs),length(ells));

for indlvlH = 1:length(lvlHs)
    lvlH = lvlHs(indlvlH);
    lvlH

    TH = refineMesh(getMesh('Square'),lvlH);
    [Th,P,P0] = refineMesh(TH,lvlh-lvlH);
    
    Sh = assembleStiffness_CR(Th,coeff);
    Dh = - assembleDivergence_CR(Th);
    
    eHlens = sqrt(sum((TH.p(TH.e(:,1),:)-TH.p(TH.e(:,2),:)).^2,2));
    ehlens = sqrt(sum((Th.p(Th.e(:,1),:)-Th.p(Th.e(:,2),:)).^2,2));
    
    % coarse interior edges
    dofH = (sum(TH.e2t) > 2-1e-10)';
    % fine interior edges
    dofh = (sum(Th.e2t) > 2-1e-10)';
    % overall dofs
    dofsys = [dofh;dofh;true(Th.nt+1,1)];
    
    inteHind = find(dofH);
    
    % compute projection operator needed for LOD
    tmp = P*sparse(TH.e,repmat(1:TH.ne,2,1)',1,TH.np,TH.ne);
    nodeindices = tmp>1-1e-10;
    p2e = sparse(repmat(1:Th.ne,2,1)',Th.e,1,Th.ne,Th.np);
    tmp2 = p2e*nodeindices;
    eHh = tmp2>2-1e-10;
    Projop = spdiags(1./eHlens,0,TH.ne,TH.ne)*(eHh'*spdiags(ehlens,0,Th.ne,Th.ne));
    Projopsys = blkdiag(Projop,Projop);
    
    % assemble global matrix for LOD saddle point problem
    Shsys = blkdiag(Sh,Sh);
    % projection operator
    Gh = spdiags(simpvol(TH)./sum(P0)',0,TH.nt,TH.nt)*P0';
    
    Asys = [Shsys                 Dh'       Projopsys'        sparse(2*Th.ne,TH.nt);
            Dh                    sparse(Th.nt,Th.nt+2*TH.ne) Gh';
            Projopsys             sparse(2*TH.ne,Th.nt+2*TH.ne+TH.nt);
            sparse(TH.nt,2*Th.ne) Gh        sparse(TH.nt,2*TH.ne+TH.nt)];
    
    %% basis functions of prototypical LOD
    % loop over all patches around interior edges
    patches = ones(TH.nt);
    lodfuns1 = sparse(2*Th.ne,length(inteHind));
    lodfuns2 = sparse(2*Th.ne,length(inteHind));
    
    inte2tH = TH.e2t(:,dofH);
    
    for indpatch = 1:length(inteHind)
        % coarse patch elements
        indtHpatch = logical(patches*inte2tH(:,indpatch));
        indthpatch = logical(P0*indtHpatch);
    
        % coarse dofs of patch
        dofHp = (sum(TH.e2t(indtHpatch,:)) > 2-1e-10)';
        % fine dofs of patch
        dofhp = (sum(Th.e2t(indthpatch,:)) > 2-1e-10)';
    
        % dofs of patch problem
        dofsys = [dofhp;dofhp;indthpatch;dofHp;dofHp;indtHpatch];
        % rhs
        tmp = zeros(TH.ne,1);
        tmp(inteHind(indpatch)) = 1;
        rhs = [zeros(2*Th.ne,2);
               zeros(Th.nt,2);
               kron(eye(2),tmp);
               zeros(TH.nt,2)];
        % solve linear system of equations
        x = zeros(length(dofsys),2);
        x(dofsys,:) = Asys(dofsys,dofsys)\rhs(dofsys,:);
        protlodfuns1(:,indpatch) = x(1:2*Th.ne,1);
        protlodfuns2(:,indpatch) = x(1:2*Th.ne,2);
    end % for
        
    protlodfuns = [protlodfuns1 protlodfuns2];
        
    %% basis functions of LOD    
    for indell = 1:length(ells)
        ell = ells(indell);
        ell
    
        % loop over all patches around interior edges
        patches = getPatches(TH,ell);
        lodfuns1 = sparse(2*Th.ne,length(inteHind));
        lodfuns2 = sparse(2*Th.ne,length(inteHind));
        
        inte2tH = TH.e2t(:,dofH);
        
        for indpatch = 1:length(inteHind)
            % coarse patch elements
            indtHpatch = logical(patches*inte2tH(:,indpatch));
            indthpatch = logical(P0*indtHpatch);
        
            % coarse dofs of patch
            dofHp = (sum(TH.e2t(indtHpatch,:)) > 2-1e-10)';
            % fine dofs of patch
            dofhp = (sum(Th.e2t(indthpatch,:)) > 2-1e-10)';
        
            % dofs of patch problem
            dofsys = [dofhp;dofhp;indthpatch;dofHp;dofHp;indtHpatch];
            % rhs
            tmp = zeros(TH.ne,1);
            tmp(inteHind(indpatch)) = 1;
            rhs = [zeros(2*Th.ne,2);
                   zeros(Th.nt,2);
                   kron(eye(2),tmp);
                   zeros(TH.nt,2)];
            % solve linear system of equations
            x = zeros(length(dofsys),2);
            x(dofsys,:) = Asys(dofsys,dofsys)\rhs(dofsys,:);
            lodfuns1(:,indpatch) = x(1:2*Th.ne,1);
            lodfuns2(:,indpatch) = x(1:2*Th.ne,2);
        end % for
        
        lodfuns = [lodfuns1 lodfuns2];
        
        %% compute errors
        e = lodfuns - protlodfuns;
        locerrH1(indlvlH,indell) = max(sqrt(diag(e'*Shsys*e)));
        save('localization.mat','locerrH1');
    end % for
end % for