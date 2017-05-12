function varargout = process_scan_subspace( varargin )
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
sProcess.Comment     = 'Scan signal subspace';
sProcess.FileTag     = '';
sProcess.Category    = 'File';
sProcess.SubGroup    = 'Sources';
sProcess.Index       = 113;
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data'};
sProcess.OutputTypes = {'results'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
sProcess.isSeparator = 0;

% === 
sProcess.options.label3.Type    = 'label';
sProcess.options.label3.Comment = '<HTML><BR>Select dimensions based on:';
sProcess.options.rsspace.Comment = {'number of dimensions','p-value threshold','SNR values'};
sProcess.options.rsspace.Type    = 'radio';
sProcess.options.rsspace.Value   = 0;
% ===
sProcess.options.ndim.Comment = 'Number of dimensions:';
sProcess.options.ndim.Type    = 'value';
sProcess.options.ndim.Value   = {4, '', 0};
% ===
sProcess.options.alpha.Comment = 'p-value threshold:';
sProcess.options.alpha.Type    = 'value';
sProcess.options.alpha.Value   = {.05, '', 3};
% ===
sProcess.options.snr.Comment = 'SNR threshold:';
sProcess.options.snr.Type    = 'value';
sProcess.options.snr.Value   = {3, '', 3};
% === Order of components
sProcess.options.label2.Type    = 'label';
sProcess.options.label2.Comment = '<HTML><BR>Order:';
sProcess.options.order.Comment = {'first', 'last'};
sProcess.options.order.Type    = 'radio';
sProcess.options.order.Value   = 0;
% === Number of scan and project iterations
sProcess.options.niter.Comment = 'Number of scanning iterations';
sProcess.options.niter.Type    = 'value';
sProcess.options.niter.Value   = {1, '', 0};
% ===
sProcess.options.iseed.Comment = 'Seed vertex for first iteration:';
sProcess.options.iseed.Type    = 'value';
sProcess.options.iseed.Value   = {0, '', 0};
end



%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFile = Run(sProcess, sInput)

% Get options
Options.islast = sProcess.options.order.Value-1;
Options.nDim = sProcess.options.ndim.Value{1};
Options.p_crit = sProcess.options.alpha.Value{1};
Options.snr_crit = sProcess.options.snr.Value{1};
Options.dimselect = sProcess.options.rsspace.Value;
Options.nIter = sProcess.options.niter.Value{1};
if isfield(sProcess.options,'issave')&&~(sProcess.options.issave.Value)
    Options.isSave = 0;
else
    Options.isSave = 1;
end
iSeed = sProcess.options.iseed.Value{1};

% Load necessary files
sInputs = bst_process('GetInputStruct', sInput.FileName);
fprintf('Subject %s\n', sInputs.SubjectName);
sMat = in_bst_data(sInput.FileName);
ChannelMat = in_bst_channel(sInputs.ChannelFile);
sStudy = bst_get('Study', sInput.iStudy);
GainFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;
GainMat = in_bst_headmodel(GainFile);
iGoodChan = sMat.GoodChannel;
NoiseCov = in_bst(sStudy.NoiseCov.FileName);

% Load subspace
if ~Options.islast
    if Options.dimselect==1
        nDim = Options.nDim;
    elseif Options.dimselect==2
        nDim = find(diff(cumsum(sMat.pval'<Options.p_crit))==0,1);
        if (nDim==1)&&(sMat.pval(1)>Options.p_crit)
            nDim = 0;
        end
    elseif Options.dimselect==3
        nDim = sum(sMat.SNR>Options.snr_crit);
    end
    Us = sMat.F(iGoodChan,1:nDim);
else
    if Options.dimselect==1
        nDim = Options.nDim;
    elseif Options.dimselect==2
        nDim = find(diff(cumsum(flip(sMat.pval)'<Options.p_crit))==0,1);
        if (nDim==1)&&(sMat.pval(1)>Options.p_crit)
            nDim = 0;
        end
    elseif Options.dimselect==3
        nDim = sum(sMat.SNR<Options.snr_crit);
    end
    Us = sMat.F(iGoodChan,end-nDim+1:end);
end

% Load leadfield
G = GainMat.Gain(iGoodChan,:);
nDip = size(G,2)/3;
G = reshape(G, size(G,1), 3, nDip);

% Apply SSP projector to leadfields
if ~isempty(ChannelMat.Projector)
    tmp = process_ssp2('BuildProjector',ChannelMat.Projector,1);
    if ~isempty(tmp)
        Proj = tmp(iGoodChan,iGoodChan);
        for iDip = 1:size(G,3)
            G(:,:,iDip) = Proj * G(:,:,iDip);
        end
    end
end

% Apply projection from subspace computation to leadfields
if isfield(sMat.Options, 'isProj')&&sMat.Options.isProj %isfield(sMat.Options, 'proj')&&~isempty(sMat.Options.proj)
    if ~isempty(ChannelMat.Projector)
        error('No!');
    end
    proj = sMat.Options.proj;
    proj = proj ./ std(proj);
    Proj = eye(length(proj))-proj*pinv(proj'*proj)*proj';
    for iDip = 1:size(G,3)
        G(:,:,iDip) = Proj * G(:,:,iDip);
    end    
end

% Apply whitener to gain and subspace
Gorig = G;
Usorig = Us;
G = [];
C = NoiseCov.NoiseCov(iGoodChan,iGoodChan);
[U, S] = svd(C,'econ');
epsilon = .001; % regularize
s = diag(S);
ss = sum(s);
sss = 1 - (cumsum(s)./ss);
i = find(sss<epsilon);
i = i(1);
Proj = sqrt(pinv(S(1:i,1:i)))*U(:,1:i)';
Us = Proj * Usorig;
for iDip = 1:nDip
    G(:,:,iDip) = Proj * Gorig(:,:,iDip);
end

% Scan
sc = scan_subspace(Us, G, Options.nIter, iSeed);

% Save source file
Options.nDim = nDim;
if ~Options.islast
    Comment = ['Subscan (plus, ' num2str(Options.nDim) ')' sMat.Comment(length('Subspace patterns')+1:end)];
else
    Comment = ['Subscan (minus, ' num2str(Options.nDim) ')' sMat.Comment(length('Subspace patterns')+1:end)];
end
ResultMat = db_template('resultsmat');
ResultMat.ImageGridAmp = sc;
ResultMat.Time = 1:size(sc,2);
ResultMat.Comment = Comment;
ResultMat.HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;
ResultMat.Options = Options;
ResultMat.Options.SubspaceFile = sInputs.FileName;
ResultMat.HeadModelType = GainMat.HeadModelType;
ResultMat.SurfaceFile = GainMat.SurfaceFile;
ResultMat.GridLoc = GainMat.GridLoc;
ResultMat.nAvg = 1;
ResultMat.Function = 'subcorr';
ResultMat.Options.iSeed = iSeed;

if Options.isSave
    OutputFile=bst_process('GetNewFilename',bst_fileparts(sInput.FileName), 'results');
    OutputFile = file_unique(OutputFile);
    bst_save(OutputFile, ResultMat);
    db_add_data(sInput.iStudy, OutputFile, ResultMat);
else
    ResultMat.Options.G = G;
    ResultMat.Options.Us = Us;
    ResultMat.Options.Gorig = Gorig;
    ResultMat.Options.Usorig = Usorig;
    ResultMat.Options.iGoodChan = iGoodChan;
    OutputFile = ResultMat;
end

end

