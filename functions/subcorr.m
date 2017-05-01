function varargout = subcorr(A, B, dim)

% if min(size(A))~=1
%     error('A must be vector');
% end

if nargin<3
    dim = 1;
end

if nargout==1
    [Ua,~]=svd(A,'econ');
    [Ub,Sb]=svd(B,'econ');
    Sc=svd(Ua'*Ub,'econ');
    if dim==1
        varargout{1} = max(Sc);       
    else
        scsort = sort(Sc,'descend');
        varargout{1} = mean(scsort(1:dim));
    end
elseif nargout>1
    [Ua, Sa, Va]=svd(A,'econ');
    [Ub, Sb, Vb]=svd(B,'econ');
    [Uc, Sc, Vc]=svd(Ua'*Ub,'econ');    
    X = Va*pinv(Sa)*Uc;
    Y = Vb*pinv(Sb)*Vc;
    if dim==1
        [varargout{1}, iMax] = max(diag(Sc));
        varargout{2} = X(:,1);
        varargout{3} = Y(:,1);
    else
        [scsort, iSort] = sort(Sc,'descend');
        varargout{2} = X(:,iSort(1,1:dim));
        varargout{3} = Y(:,iSort(1,1:dim));
    end
end

end


