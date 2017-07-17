function [p_maj, gamma_0, p_kk_gamma, gammavals, p_k_gamma] = prevalence_test(observed, alpha_ind, beta_ind, alpha_group)

% [p_maj, gamma_0, p_kk_gamma, gammavals, p_k_gamma] = prevalence_test(input, alpha_ind, beta_ind, alpha_group)
%
% Computes prevalence test
%
% Inputs
%
% observed: either vector of [n_obs, n_total] (n_obs: number of subjects with effect, 
%                                           n_total: number of subjects),
%        or vector of p-values [p_subj1, p_subj2, ..., p_subjN]
% alpha_ind: p-value threshold, if observed is vector of p-values
% beta_ind: sensitivity (probably set to 1)
% alpha_group: see Outputs
%
% Outputs
%
% p_maj: p-value for majority hypothesis
% gamma_0: highest gamma_0 that can be rejected at threshold of alpha_group

% Is input p-values or k out of n?
if sum(observed>1)==0
    pvals = observed;
else
    pvals = [zeros(observed(1),1); ones(observed(2)-observed(1),1)];
end
if nargin<2
    alpha_ind = .05;
end
if nargin<3
    beta_ind = 1;
end
if nargin<4
    alpha_group = .05;
end

dgamma = .001;

n = size(pvals,1);
k_obs = sum(pvals<alpha_ind);
gammavals = 0:dgamma:1;
kvals = 0:n;

clear p_k_gamma tmp
for iGam = 1:numel(gammavals)
    p_pos = gammavals(iGam).*beta_ind + (1-gammavals(iGam)).*alpha_ind;
    p_neg = gammavals(iGam).*(1-beta_ind) + (1-gammavals(iGam)).*(1-alpha_ind);
    for iK = 1:numel(kvals)
        k = kvals(iK);
        tmp(iK) = nchoosek(n,k)*p_pos^k*p_neg^(n-k);
    end
    p_k_gamma(iGam,:) = tmp; % k or more subjects out of n
    p_kk_gamma(iGam,:) = 1-cumsum(tmp)+tmp; % k or more subjects out of n
end

iK = find(kvals==k_obs);
p_maj = p_kk_gamma(gammavals==.5,iK);

index = sum(p_kk_gamma(:,iK) < alpha_group);
if index~=0
    gamma_0 = gammavals(index);
else
    gamma_0 = 0;
end

end


