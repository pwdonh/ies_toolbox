function [SNR, Results, C1, C2, C1_all, C2_all, Pow_all, Pow_all_unreg, C_orig, y_all] = ...
    compute_signal_subspace_core(X1, X2, y, pipeline, Options)  

% Compute covariance matrices
% fprintf('Compute covariance\n');
if strcmp(pipeline, 'pca') 
    % if PCA, C2 reduces to the identity
    [C1, C1_all] = compute_covariance(X1, 1);
    C2 = Options.NoiseCov;
    C2_all = repmat(C2,1,1,size(C1_all,3));
elseif strcmp(pipeline, 'coh_ref')
    % if Coherence, the signal and noise spaces are defined using the
    % cross- and auto-correlation functions
    error('Not implemented.')
    [C1, C2, C1_all, C2_all] = compute_auto_cross_covariance(X1, y, Options.PassBand, Options.fs);
elseif strcmp(pipeline, 'coh_ref_am')
    [C1_all, C2_all, y_all, C_orig] = compute_cov_am(X1, y, 1, Options.nBins);
    C1 = mean(C1_all,3);
    C2 = mean(C2_all,3);
elseif strcmp(pipeline, 'aec')
    [C1_all, C2_all, y_all, C_orig] = compute_cov_am(X1, y, Options.TypeAM, Options.nBins);
    C1 = mean(C1_all,3);
    C2 = mean(C2_all,3);
else % Condition 1 over 2, Frequency band 1 over 2, etc.
    [C1, C1_all] = compute_covariance(X1, 1);
    [C2, C2_all] = compute_covariance(X2, 1);
end

if ~strcmp(pipeline, 'coh_ref_am')&&~strcmp(pipeline, 'aec')
    C_orig = [];
    y_all = [];
end

% Project out given topography
if Options.isProj %isfield(Options, 'proj')&&~isempty(Options.proj)
    proj = Options.proj ./ std(Options.proj);
    P = eye(length(proj))-proj*pinv(proj'*proj)*proj';
    C1 = P*C1*P'; C2 = P*C2*P';
    for iSegment = 1:size(C1_all,3)
        C1_all(:,:,iSegment) = P*C1_all(:,:,iSegment)*P';
        C2_all(:,:,iSegment) = P*C2_all(:,:,iSegment)*P';
    end
end

% Regularize covariance matrices
C1_all_unreg = C1_all; C2_all_unreg = C2_all;
[C1, C1_all] = regularize_covariance(C1, Options.alpha, C1_all);
[C2, C2_all] = regularize_covariance(C2, Options.alpha, C2_all);

% Solve generalized eigenvalue problem
% fprintf('Compute subspace\n');
for iSegment = 1:size(C1_all,3)
    inputs1{iSegment} = C1_all(:,:,iSegment);
    inputs2{iSegment} = C2_all(:,:,iSegment);
end  
[SNR, Results{2}, Results{1}] = geneig_func({inputs1, inputs2}, Options);
nTrial = size(C1_all,3);
for iTrial = 1:nTrial
    Pow_all(:,iTrial,1) = diag(Results{2}'*C1_all(:,:,iTrial)*Results{2});
    Pow_all(:,iTrial,2) = diag(Results{2}'*C2_all(:,:,iTrial)*Results{2});
    Pow_all_unreg(:,iTrial,1) = diag(Results{2}'*C1_all_unreg(:,:,iTrial)*Results{2});
    Pow_all_unreg(:,iTrial,2) = diag(Results{2}'*C2_all_unreg(:,:,iTrial)*Results{2});    
end

end