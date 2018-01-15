function varargout = process_image_subspace( varargin )
%
%
%
%

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
% Authors: Peter Donhauser 2015

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
% Description the process
sProcess.Comment     = 'Image signal subspace';
sProcess.FileTag     = '';
sProcess.Category    = 'File';
sProcess.SubGroup    = 'iES plugin';
sProcess.Index       = 113;
% Definition of the input accepted by this process
sProcess.InputTypes  = {'results'};
sProcess.OutputTypes = {'results'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
sProcess.isSeparator = 0;

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFile = Run(sProcess, sInput)

% Load imaging kernel
KernelMat = in_bst_results(sInput.FileName);
K = KernelMat.ImagingKernel;

% Load covariance matrices
sMat = in_bst_data(KernelMat.DataFile);

% Compute metric of interest at each source location
Kr = reshape(K', size(K,2), 3, size(K,1)/3);
for iSource = 1:size(Kr,3)
    C1 = Kr(:,:,iSource)'*sMat.C1*Kr(:,:,iSource);
    C2 = Kr(:,:,iSource)'*sMat.C2*Kr(:,:,iSource);
    if strcmp(sMat.Pipeline,'aec')
        [U,S]=svd(C2);
        [~,ii]=max(abs(diag(S)));
        Ko(iSource,:) = Kr(:,:,iSource)*U(:,ii);
    else
%         geneigs = geneig_func({{C1},{C2}});
%         [~,ii]=max(abs(geneigs));
%         y(iSource,1) = geneigs(ii);
        y(iSource,1) = max(svd(C1))/max(svd(C2));
    end
end
if strcmp(sMat.Pipeline,'aec')
    y = compute_correlation_from_cov(sMat.C1, sMat.C2_all, Ko')';
end

% Prepare result structure
ResultMat = KernelMat;
ResultMat.ImagingKernel = [];
ResultMat.ImageGridAmp = [y y];
ResultMat.Time = 1:2;
ResultMat.Comment = 'Subspace imaged';
ResultMat.nComponents = 1;
ResultMat.Options.SubspaceFile = KernelMat.DataFile;

% save file
OutputFile = bst_process('GetNewFilename', bst_fileparts(KernelMat.DataFile), 'results');
OutputFile = file_unique(OutputFile);
bst_save(OutputFile, ResultMat);
db_add_data(sInput.iStudy, OutputFile, ResultMat);

end
