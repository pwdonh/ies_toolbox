function [C_r, C_all_r] = regularize_covariance(C, alpha, C_all)

% Add to diagonal without changing sum(svd(C))
n = size(C,2);
C_r = (1 - alpha) * C; 
mu = trace(C)/n; 
C_r = C_r + eye(n)*mu*alpha;
if nargin>2
    % Also for single segment C's, so that average C is the same
    C_all_r = bsxfun(@plus,(1-alpha)*C_all,eye(n)*mu*alpha);
else
    C_all_r = [];
end

end