function [e,te,e2t,p2e] = getEdges(mesh)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get edges and information about their connectivity
    %
    % Input:
    %     mesh:  simplicial mesh 
    %   
    % Output: 
    %        e:  tuples of node indices of edges
    %       te:  triples of edge indices of elements
    %      e2t:  connectivity matrix between elements and edges
    %      p2e:  connectivity amtrix between nodes and edges
    %
    % M. Hauck, A. Lozinski
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$

    % extracts sets of edges 
    t = mesh.t;
    e1 = t(:,[2 3]);
    e2 = t(:,[3 1]);
    e3 = t(:,[1 2]);

    % as sets of their nodes
    vs = zeros(size(t,1)*3,2);
    vs(1:3:end,:) = e1;
    vs(2:3:end,:) = e2;
    vs(3:3:end,:) = e3;

    % repeated sets of nodes (joint edges) are eliminated 
    [e,te] = deleterepeatedrows(vs);
    te = reshape(te,size(t,2),size(t,1))';

    % incidence matrices
    if nargout > 2
        e2t = sparse(repmat(1:size(te,1),3,1)',te,1,size(te,1),size(e,1));
    end % if
    if nargout > 3
        p2e = sparse(repmat(1:size(e,1),2,1)',e,1,size(e,1),max(max(e)));
    end % if
end % function

function r = insertvector(v,pos_from,pos_to)
     tf = false(1,numel(v)+numel(pos_to));
     r = double(tf);
     tf(pos_to) = true;
     r(~tf) = v;
     r(tf) = r(pos_from);    
end % function

function [matrix,I] = deleterepeatedrows(matrix)
    [matrixs,tags] = sortrows(sort(matrix,2));
    k = find(all(diff(matrixs)==0,2));
    % these rows of matrix are repeated 
    repeated = tags(k);
    % and these rows will be removed
    removed = tags(k+1);
    % both lists are sorted 
    [removeds, tags2] = sort(removed);
    repeateds = repeated(tags2);
    % delete the tags to removed rows
    tags(k+1) = [];
    % and recover the original array, in the original order.
    matrix = matrix(sort(tags),:);
    % row indices before matrix compression indicating repetition   
    I = insertvector((1:size(matrix,1))',repeateds,removeds);
end % function