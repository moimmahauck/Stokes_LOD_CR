%% This code reproduces the second numerical experiment in the paper
%
% A Localized Orthogonal Decomposition Method for Heterogeneous Stokes
% Problems
%
% Note: the below parameter configuration was used on a cluster. The
% computation may not be feasuble on a laptop.
%
% M. Hauck, A. Lozinski

init;
clear all; close all; 

parpool("Threads")

ells = 2:7;
lvlHs = 1:6;

lvlh = 10;

% rhs
f = @(x) [x(:,2) -x(:,1)];

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
coeff = P0*coeff;

%% compute reference solution
dofh = (sum(Th.e2t) > 2-1e-10)';
Sh = assembleStiffness_CR(Th,coeff);
Dh = - assembleDivergence_CR(Th);
Mh = assembleMass_CR(Th);
Ch = simpvol(Th);
P = computeCGtoCR(Th);

% incorporating zero mean constraint to make system aolvable
dofsys = [dofh;dofh;true(Th.nt+1,1)];
Asys = [blkdiag(Sh,Sh)    Dh'                 sparse(2*Th.ne,1);
        Dh                sparse(Th.nt,Th.nt) Ch;
        sparse(1,2*Th.ne) Ch'                 sparse(1,1)];
rhs1 = reshape(Mh*P*reshape(f(Th.p),[],2),[],1);
rhs2 = zeros(Th.nt,1);
rhssys = [rhs1;rhs2; 0];
x = zeros(length(dofsys),1);
x(dofsys) = Asys(dofsys,dofsys)\rhssys(dofsys);
uh = x(1:2*Th.ne);
uh1 = x(1:Th.ne); 
uh2 = x(Th.ne+1:2*Th.ne);
ph = x(2*Th.ne+1:end-1);

%% basis functions of LOD
errL2ref_u = zeros(length(ells),length(lvlHs));
errH1ref_u = zeros(length(ells),length(lvlHs));
errL2ref_p = zeros(length(ells),length(lvlHs));
errL2ref_pp = zeros(length(ells),length(lvlHs));

fprintf('rough: lvlHmin = %d, lvlHmax = %d, ellmin = %d, ellmax = %d, lvlh = %d.\n',lvlHs(1),lvlHs(end),ells(1),ells(end),lvlh);
for indlvlH = 1:length(lvlHs)
    lvlH = lvlHs(indlvlH);
    lvlH        

    global_flag = false;
    for indell = 1:length(ells)
        ell = ells(indell);
        ell

        % skip iteration if already computed
        if global_flag
            fprintf('skip.\n');
            errL2ref_u(indell,indlvlH) = errL2ref_u(indell-1,indlvlH);
            errH1ref_u(indell,indlvlH) = errH1ref_u(indell-1,indlvlH)
            errL2ref_p(indell,indlvlH) = errL2ref_p(indell-1,indlvlH)
            errL2ref_pp(indell,indlvlH) = errL2ref_pp(indell-1,indlvlH)
            continue;
        end % if

        TH = refineMesh(getMesh('Square'),lvlH);
        [Th,P,P0] = refineMesh(TH,lvlh-lvlH);
        eHlens = sqrt(sum((TH.p(TH.e(:,1),:)-TH.p(TH.e(:,2),:)).^2,2));
        ehlens = sqrt(sum((Th.p(Th.e(:,1),:)-Th.p(Th.e(:,2),:)).^2,2));
        
        % coarse interior edges
        dofH = (sum(TH.e2t) > 2-1e-10)';
        % fine interior edges
        dofh = (sum(Th.e2t) > 2-1e-10)';
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
        Mhsys = blkdiag(Mh,Mh);
        Shsys = blkdiag(Sh,Sh);
        % projection operator
        Gh = spdiags(simpvol(TH)./sum(P0)',0,TH.nt,TH.nt)*P0';
        
        Asys = [Shsys                 Dh'       Projopsys'        sparse(2*Th.ne,TH.nt);
                Dh                    sparse(Th.nt,Th.nt+2*TH.ne) Gh';
                Projopsys             sparse(2*TH.ne,Th.nt+2*TH.ne+TH.nt);
                sparse(TH.nt,2*Th.ne) Gh        sparse(TH.nt,2*TH.ne+TH.nt)];
        
        
        % loop over all patches around interior edges
        patches = getPatches(TH,ell);
        inte2tH = TH.e2t(:,dofH);

        % detect if all patches are global
        tmp = logical(patches*inte2tH);
        global_flag = sum(logical(tmp),'all')-numel(tmp) == 0;
        
        lodfuns1 = sparse(2*Th.ne,length(inteHind));
        lodfuns2 = sparse(2*Th.ne,length(inteHind));
        lodpfuns1 = sparse(Th.nt,length(inteHind));
        lodpfuns2 = sparse(Th.nt,length(inteHind));
        
        parfor indpatch = 1:length(inteHind)
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
            lodpfuns1(:,indpatch) = x(2*Th.ne+1:2*Th.ne+Th.nt,1);
            lodpfuns2(:,indpatch) = x(2*Th.ne+1:2*Th.ne+Th.nt,2);
        end % for
        
        lodfuns = [lodfuns1 lodfuns2];
        lodpfuns = [lodpfuns1 lodpfuns2];

        %% LOD method
        Slod = lodfuns'*(Shsys*lodfuns);
        Dlod = P0'*Dh*lodfuns;
        Clod = P0'*Ch;
        
        Alod = [Slod                   Dlod'               sparse(size(Slod,1),1);
                Dlod                   sparse(TH.nt,TH.nt) Clod;
                sparse(1,size(Slod,1)) Clod'               sparse(1,1)];
        
        rhslod = [lodfuns'*rhs1; P0'*rhs2; 0];
        
        xlod = Alod\rhslod;
        ulod = lodfuns*xlod(1:sum([dofH;dofH]));
        ulod1 = ulod(1:Th.ne);
        ulod2 = ulod(Th.ne+1:end);
        plodP0 = xlod(end-TH.nt:end-1);
        plodpp = P0*plodP0 + lodpfuns*xlod(1:sum([dofH;dofH]));

        %% checks for lod solution
        e_u = uh-ulod;
        e_p = spdiags(1./sum(P0)',0,TH.nt,TH.nt)*P0'*ph-plodP0;
        e_pp = ph - plodpp;
                
        % H1-error
        errL2ref_u(indell,indlvlH) = sqrt(e_u'*(Mhsys*e_u));
        errH1ref_u(indell,indlvlH) = sqrt(e_u'*(Shsys*e_u))
        errL2ref_p(indell,indlvlH) = sqrt(e_p'*(simpvol(TH).*e_p))
        errL2ref_pp(indell,indlvlH) = sqrt(e_pp'*(simpvol(Th).*e_pp))
        save(['convergence' num2str(lvlh) '.mat'],'errL2ref_u','errH1ref_u','errL2ref_p','errL2ref_pp','ells','lvlHs','lvlh');
    end % for
end % for