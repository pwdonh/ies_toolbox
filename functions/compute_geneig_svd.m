function varargout = compute_geneig_svd(CovA, CovB, epsilon)

% [SNR, Filters, Patterns] = compute_geneig_svd(CovA, CovB, epsilon)
%
% Output: SNR (generalized eigenvalues)
%         Filters
%         Patterns

[UB, SB, ~] = svd(CovB,'econ');
s = diag(SB);
ss = sum(s);
sss = 1 - (cumsum(s)./ss);
i = find(sss<epsilon);
if ~isempty(i)
    i = i(1);
else
    i = numel(s);
end
si = 1./s;
Proj = diag(sqrt(si(1:i)))*UB(:,1:i)';
CovA_w = Proj*CovA*Proj';
[UA_w, SA_w, ~] = svd(CovA_w,'econ');

% Generalized eigenvalues
varargout{1} = diag(SA_w);
if nargout>1
    % Reverse whitening operation
    varargout{2} = columnnorm(UB(:,1:i)*diag(sqrt(si(1:i)))*UA_w);    
    varargout{3} = columnnorm(UB(:,1:i)*diag(sqrt(s(1:i)))*UA_w);
    [varargout{1}, index] = sort(diag(varargout{2}'*CovA*varargout{2}) ./ ...
        diag(varargout{2}'*CovB*varargout{2}),'descend');
    varargout{2} = varargout{2}(:,index);
    varargout{3} = varargout{3}(:,index);
end

end









% function varargout = compute_geneig_matlab(CovA, CovB)
% 
% % Solve generalized eigenvalue problem
% [V, D] = eig(CovA, CovB);
% D = diag(D);
% [~, index] = sort(D, 'descend');
% SNR = D(index);
% varargout{2} = V(:, index);
% % invert the filters to obtain the forward patterns
% if nargout>2
%     varargout{3} = pinv(varargout{2})';
% end
% % Other output
% varargout{1} = SNR;
% 
% end