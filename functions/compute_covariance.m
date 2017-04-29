function [C, C_all] = compute_covariance(X, isKeepTrials)

if nargin<2
    isKeepTrials = 0;
end

X = bsxfun(@minus, X, mean(X,2));
nTrial = size(X,3);
nSample = size(X,2)*nTrial;
denom = nSample-1;

if ~isKeepTrials
    % Compute by accumulating matrix product trial by trial and dividing
    % by the number of samples
    C_all = [];
    C = zeros(size(X,1));
    for iTrial = 1:size(X,3)
        C = C + X(:,:,iTrial) * X(:,:,iTrial)';
    end
    C = C ./ denom;
else
    % Keep all
    C_all = zeros(size(X,1),size(X,1),size(X,3));
    for iTrial = 1:size(X,3)
        C_all(:,:,iTrial) = X(:,:,iTrial) * X(:,:,iTrial)';
    end
    C_all = C_all ./ (denom/nTrial);
    C = sum(C_all,3) ./ nTrial;    
end

end