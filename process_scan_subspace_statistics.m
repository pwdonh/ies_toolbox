function varargout = process_scan_subspace_statistics( varargin )
% 
%
% USAGE:  

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
%
% Copyright (c)2000-2014 University of Southern California & McGill University
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
% Authors: Peter Donhauser

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
% Description the process
sProcess.Comment     = 'Subspace scanning group statistics';
sProcess.FileTag     = '';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Sources';
sProcess.Index       = 302;
% Definition of the input accepted by this process
sProcess.InputTypes  = {'results'};
sProcess.OutputTypes = {'results'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
sProcess.isSeparator = 1;

% 
sProcess.options.nperm.Comment = 'Number of permutation:';
sProcess.options.nperm.Type    = 'value';
sProcess.options.nperm.Value   = {1000, '', 0};
% 
sProcess.options.isinterp.Comment = 'Run permutation test in MRI volume';
sProcess.options.isinterp.Type    = 'checkbox';
sProcess.options.isinterp.Value   = 0;
%
sProcess.options.niter.Comment = 'Number of scanning iterations:';
sProcess.options.niter.Type    = 'value';
sProcess.options.niter.Value   = {1, '', 0};

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput)

nPerm = sProcess.options.nperm.Value{1};
isInterp = sProcess.options.isinterp.Value;
nIter = sProcess.options.niter.Value{1};

sProcessAvg = process_average('GetDescription');
sProcessAvg.options.avgtype.Value = 1; % Everything
AvgFile = process_average('Run', sProcessAvg, sInput);

sProcessScan = process_scan_subspace('GetDescription');
sProcessScan.options.rsspace.Value = 1; % Set dimensions manually
sProcessScan.options.niter.Value{1} = nIter;
sProcessScan.options.issave.Value = 0;

%%

if isInterp
    sMat = in_bst_results(sInput(1).FileName);
    sSubject = bst_get('Subject', sInput(1).SubjectFile);
    sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);
    grid2mri_interp = grid_interp_mri(sMat.GridLoc, sMri, sMat.SurfaceFile, 0);
end
    
%%

for iFile = 1:numel(sInput)
    
    fprintf('Loading map %d, ', iFile);
    sMat = in_bst_results(sInput(iFile).FileName);
    sProcessScan.options.ndim.Value{1} = sMat.Options.nDim;
    if isfield(sMat.Options, 'iSeed')
        sProcessScan.options.iseed.Value{1} = sMat.Options.iSeed;
    end
    if sMat.Options.islast
        sProcessScan.options.order.Value = 1; % first
    else
        sProcessScan.options.order.Value = 2; % first
    end
    sInputScan = bst_process('GetInputStruct', sMat.Options.SubspaceFile);
    sMatReversed = process_scan_subspace('Run', sProcessScan, sInputScan);
    data = sMat.ImageGridAmp(:,nIter);
    data_reversed = sMatReversed.ImageGridAmp(:,nIter);
    
    if isInterp
        data = tess_interp_mri_data(grid2mri_interp, size(sMri.Cube), data, 1);
        data_reversed = tess_interp_mri_data(grid2mri_interp, size(sMri.Cube), data_reversed, 1);
    end
        
    % save subcorr maps
    scMaps{1}{iFile} = data(:);
    scMaps{2}{iFile} = data_reversed(:);
    
end

p_crit = [.05 .01 .005 .001];

% Perform permutation test
Options.nPerm = nPerm;
Options.nDim = 1;
Options.stat0 = 0;
Options.alpha = .05;
Options.verbose = 1;
Options.func = @(inputs,para) mean([inputs{1}{:}],2);
Options.isMaximum = 1;
% Options.two_tailed = 0; % one-tailed since only positive values
[pval_max, stat_obs, permvals] = permtest_dep2_new(scMaps{1}, scMaps{2}, Options);
% sortperm = sort(max(permvals));
sortperm = sort(permvals);
for iCrit = 1:numel(p_crit)
    ii = ceil((1-(p_crit(iCrit)/2))*size(permvals,2));
    stat_crit(iCrit) = sortperm(ii);
end

Results.stat_obs = stat_obs;
Results.pval_max = pval_max;
Results.max_permvals = max(permvals);
Results.stat_crit = stat_crit;
Results.p_crit = p_crit;
Results.InputFiles = {sInput.FileName};

AvgMat = in_bst_results(AvgFile{1});
AvgMat.ImageGridAmp = [AvgMat.ImageGridAmp(:,nIter) AvgMat.ImageGridAmp(:,nIter)];
AvgMat.StatResults = Results;
AvgMat.StatResults.ptype = 'pval_max';
AvgMat.StatResults.isinterp = isInterp;
bst_save(AvgFile{1}, AvgMat);

OutputFiles = AvgFile;

end
























