function varargout = geneig_func(inputs, para)

% Output: SNR (generalized eigenvalues)
%         Filters
%         Patterns

if nargin<2
    para.epsilon = 0;
end
if ~isfield(para,'epsilon')
    para.epsilon = 0;
end

CovA = mean(cat(3,inputs{1}{:}), 3);
CovB = mean(cat(3,inputs{2}{:}), 3);

if nargout>1
    [varargout{1}, varargout{2}, varargout{3}] = compute_geneig_svd(CovA, CovB, para.epsilon);
else
    varargout{1} = compute_geneig_svd(CovA, CovB, para.epsilon);
end

end