function varargout = process_subspace_statistics( varargin )
% 
%
% 

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2015 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Peter Donhauser, 2015

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<*DEFNU>
    % Description the process
    sProcess.Comment     = 'Compute Signal Subspace: Statistics';
    sProcess.FileTag     = '| subsp';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 112;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
%     sProcess.isSeparator = 1;
    % Number of iterations
    sProcess.options.niter.Comment = 'Number of iterations: ';
    sProcess.options.niter.Type    = 'value';
    sProcess.options.niter.Value   = {600, '', 0};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs)

nFile = numel(sInputs);
nIter = sProcess.options.niter.Value{1};

% Check that subspace files are compatible
for iFile = 1:nFile
    sMat = in_bst_data(sInputs(iFile).FileName);
    if iFile ==1
%         warning('Uncomment indSegment part again')
        indSegment = sMat.Options.indSegment;
        FilesA = sMat.FilesA;
        pipeline = sMat.Pipeline;
    else
        if ~sum(indSegment(:) == sMat.Options.indSegment(:))
            error('Files must have same windowing');
        else
            indSegment = sMat.Options.indSegment;
        end
        for iFilesA = 1:numel(sMat.FilesA)
            if sum((FilesA{iFilesA}==sMat.FilesA{iFilesA})==0)>0
                error('FilesA must be the same');
            end
            if ~strcmp(sMat.Pipeline,pipeline)
                error('Pipelines must be the same');
            end
        end
    end
end

% Compute bootstrap statistics over regularization threshold
pval_reg = {}; ci_reg = {};
for iFile = 1:nFile
    sMat = in_bst_data(sInputs(iFile).FileName);
    if sMat.Options.alpha>0
        [pval_reg{iFile}, ci_reg{iFile}] = BootstrapStatsReg(sMat, 5999);
    else
        pval_reg{iFile} = zeros(size(sMat.SNR));
        ci_reg{iFile} = [];
    end
end

% Compute statistics
switch pipeline
    case 'stim_base'
        Results = PermutationStats({sInputs.FileName}, nIter);
        pval = {Results.pval_max};
    case 'freq1_freq2'
        Results = BootstrapStats({sInputs.FileName}, nIter);
        pval = {Results.pval_snr};
    case 'coh_ref'
        Results = ShufflingStats({sInputs.FileName}, nIter);
        pval = {Results.pval_shuff};
    case 'coh_ref_am'
        Results = ShufflingStats({sInputs.FileName}, nIter);
        pval = {Results.pval_shuff};
    case 'aec'
        Results = ShufflingStats({sInputs.FileName}, nIter);
        pval = {Results.pval_shuff};        
end

% Save in file
for iFile = 1:nFile
    Results(iFile).pval_reg = pval_reg{iFile};
    Results(iFile).ci_reg = ci_reg{iFile};
    sMat = in_bst_data(sInputs(iFile).FileName);
    sMat.pval = pval{iFile};
    sMat.pval(Results(iFile).pval_reg>0) = 1;
    sMat.Results = Results(iFile);
    fullFile = file_fullpath(sInputs(iFile).FileName);
    bst_save(fullFile, sMat);
end

OutputFiles = {sInputs.FileName};

end

function [pval_reg, ci_reg] = BootstrapStatsReg(sMat, nBoot)

% Bootstrap confidence interval over regularization threshold
nSeg = size(sMat.Pow_all,2);
nChan = size(sMat.Pow_all_unreg,1);
pow1_boot = zeros(nChan, nBoot);
pow2_boot = zeros(nChan, nBoot);
for iBoot = 1:nBoot
    bootsample = randsample(1:nSeg,nSeg,1);
    pow1_boot(:,iBoot) = mean(sMat.Pow_all_unreg(:,bootsample,1),2);
    pow2_boot(:,iBoot) = mean(sMat.Pow_all_unreg(:,bootsample,2),2);
end
thresh1 = trace(sMat.C1)/nChan*sMat.Options.alpha;
thresh2 = trace(sMat.C2)/nChan*sMat.Options.alpha;
pval1 = 1-sum(pow1_boot>thresh1,2)'/nBoot;
pval2 = 1-sum(pow2_boot>thresh2,2)'/nBoot;
pow1_boot_sort = sort(pow1_boot,2);
pow2_boot_sort = sort(pow2_boot,2);
ci_reg = cat(3,pow1_boot_sort(:,[1 end]),pow2_boot_sort(:,[1 end]));
% If coh_ref_am, then first C matrix can have negative correlations
if any(any(ci_reg(:,:,1)<0))
    pval_reg = pval2;
else
    pval_reg = sum([pval1; pval2])>0;
end

end

function Results = PermutationStats(SubspaceFiles, nPerm)

% Prepare input to permutation test function
sMat = in_bst_data(SubspaceFiles{1});
nTrial = size(sMat.Pow_all,2);
C1_all = {}; C2_all = {};
for iTrial = 1:nTrial
    C1_all{iTrial} = sMat.C1_all(:,:,iTrial);
    C2_all{iTrial} = sMat.C2_all(:,:,iTrial);
end

% Permutation test based on maximum statistic
Options.nPerm = nPerm;
Options.nDim = 1;
Options.stat0 = 1;
Options.alpha = .05;
Options.verbose = 1;
Options.isMaximum = 1;
Options.func = 'geneig_func';
[pval_max, ~, stat_perm] = permtest_dep2_new(C1_all, C2_all, Options);
maxstat = max(stat_perm);

Results.pval_max = pval_max;

end

function Results = BootstrapStats(SubspaceFiles, nBoot)

% nBoot = 6000;
p_crit = .001;
nFile = numel(SubspaceFiles);
indCi = [floor(nBoot*p_crit) ceil(nBoot*(1-p_crit))];

% Load data files
for iFile = 1:nFile
    sMat = in_bst_data(SubspaceFiles{iFile});
    C1(:,:,iFile) = sMat.C1;
    C2(:,:,iFile) = sMat.C2;    
end
snr0_all = trace(sum(C1,3))./trace(sum(C2,3));

% Bootstrap confidence interval over average SNR
for iFile = 1:nFile
    sMat = in_bst_data(SubspaceFiles{iFile});
    Pow_all = sMat.Pow_all;
    nSeg = size(sMat.Pow_all,2);
    nChan = size(sMat.Pow_all,1);
    snr0 = trace(sMat.C1)/trace(sMat.C2);
    snr_boot = zeros(nChan, nBoot);
    for iBoot = 1:nBoot
        bootsample = randsample(1:nSeg,nSeg,1);
        pow1 = mean(sMat.Pow_all(:,bootsample,1),2);
        pow2 = mean(sMat.Pow_all(:,bootsample,2),2);
        snr_boot(:,iBoot) = pow1./pow2;
    end
    Results(iFile).pval_snr = 1-sum(snr_boot>snr0,2)'/nBoot;
    Results(iFile).pval_snr_all = 1-sum(snr_boot>snr0_all,2)'/nBoot;
    snr_boot_sort = sort(snr_boot,2);
    Results(iFile).ci_snr = snr_boot_sort(:,indCi);
    Results(iFile).snr0 = snr0;
    Results(iFile).snr0_all = snr0_all;
    Results(iFile).p_crit = p_crit;
end

end

function Results = ShufflingStats(SubspaceFiles, nShuff)
        
% nShuff = 100;

nFile = numel(SubspaceFiles);

% Load data files A
sMat = in_bst_data(SubspaceFiles{1});
FilesA = sMat.FilesA;
pipeline = sMat.Pipeline;
Options = sMat.Options;
Options.isSave = 0;
[~,X,~,y] = compute_signal_subspace(FilesA, sMat.FilesB, Options);
SNR_obs{1} = sMat.SNR;

% Load data files B
for iFile = 2:nFile
    sMat = in_bst_data(SubspaceFiles{iFile});
    [~,~,~,y(iFile,:,:)] = compute_signal_subspace(FilesA, sMat.FilesB, Options);
    SNR_obs{iFile} = sMat.SNR;
end

% Shuffling iterations
indSegment = Options.indSegment;
indTrial = unique(indSegment(:,2));
nTrial = length(indTrial);
nWin = size(indSegment,1);

for iShuff = 1:nShuff
    tic
    % Random shuffling indices
    shuffvec = randsample(1:nWin,nWin);
    % Compute for each file B
    for iFile = 1:nFile
        SNR = compute_signal_subspace_core(X, [], y(iFile,:,shuffvec), ...
            pipeline, Options);
        snr_max_shuff(iFile,iShuff) = SNR(1);
        snr_min_shuff(iFile,iShuff) = SNR(end);
    end
    fprintf('%d, %.2f seconds\n',iShuff, toc);
end

% Compute p-values
for iFile = 1:nFile
    pval = (sum(bsxfun(@lt, SNR_obs{iFile}', snr_max_shuff(iFile,:)),2)+1)/(nShuff+1);
    if strcmp(pipeline,'aec')||strcmp(pipeline,'coh_ref_am')
        pval_minus = (sum(bsxfun(@gt, SNR_obs{iFile}', snr_min_shuff(iFile,:)),2)+1)/(nShuff+1);
        pval(SNR_obs{1}<0) = pval_minus(SNR_obs{1}<0);
    end
    fprintf('File %d: %d significant dimensions out of %d\n', iFile, ...
        sum(pval<.05), length(pval));
    Results(iFile).pval_shuff = pval;
end

end


