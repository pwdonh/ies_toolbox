function varargout = process_compute_coregistration( varargin )
% PROCESS_MEGREG: Co-register different datasets (runs/subjects/conditions) to the same channel file.
%
% DESCRIPTION:
%     MEG runs acquired at different moments typically have different sensor positions, 
%     they cannot be averaged or compared with each other.
%     This process computes one common head position and interpolates the magnetic fields of all the
%     input files, so that they all share the same channel file at the end.
%
% USAGE:            OutputFiles = process_megreg('Run', sProcess, sInputs)
%                       maxDist = process_megreg('Compute', sInputs, isShareChan, isAvgChan, epsilon, isDebug)
%       [maxDist,ErrAvg,ErrStd] = process_megreg('Test', epsilon)


% @=============================================================================
% This function is part of the Brainstorm software:
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
% Authors: Peter Donhauser 2015
%          Sylvain Baillet 2002-2006
%          Francois Tadel, 2012

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Compute co-registration projectors';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'iES plugin';
    sProcess.Index       = 1000;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 2;
    sProcess.isSeparator = 0;
    % Definition of the options
    % === TARGET CHANNEL FILE
    sProcess.options.label1.Type    = 'label';
    sProcess.options.label1.Comment = 'Target sensors positions :';
    sProcess.options.targetchan.Comment = {'Average of all the runs', 'First channel file in the list'};
    sProcess.options.targetchan.Type    = 'radio';
    sProcess.options.targetchan.Value   = 1;
%     % === SHARE CHANNEL FILE
%     sProcess.options.label2.Type    = 'label';
%     sProcess.options.label2.Comment = '<BR>Use default channel file :';
%     sProcess.options.sharechan.Comment = {'Yes, share the same channel file between runs', 'No, do not modify the database organization'};
%     sProcess.options.sharechan.Type    = 'radio';
%     sProcess.options.sharechan.Value   = 1;
    % === EPSILON 
    sProcess.options.label3.Type    = 'label';
    sProcess.options.label3.Comment = ' ';
    sProcess.options.epsilon.Comment = 'Smoothing parameter: ';
    sProcess.options.epsilon.Type    = 'value';
    sProcess.options.epsilon.Value   = {.0001, '', 6};
    % === Reference source space
    sProcess.options.label4.Type    = 'label';
    sProcess.options.label4.Comment = '<HTML><BR>Use reference source space :';
    sProcess.options.rsspace.Comment = {'Sphere', 'Anatomical','Identity'};
    sProcess.options.rsspace.Type    = 'radio';
    sProcess.options.rsspace.Value   = 0;    
    % === Apply SSPs
    sProcess.options.applyssp.Comment = 'Apply SSPs';
    sProcess.options.applyssp.Type    = 'checkbox';
    sProcess.options.applyssp.Value   = 0;
    % === Comment
    sProcess.options.comment.Comment = 'New condition name (default @intra):';
    sProcess.options.comment.Type    = 'text';
    sProcess.options.comment.Value   = '@intra';        
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function [OutputFiles, maxDist] = Run(sProcess, sInputs) %#ok<DEFNU>

% Options
% isShareChan = (sProcess.options.sharechan.Value == 1);
isAvgChan   = (sProcess.options.targetchan.Value == 1);
epsilon     = sProcess.options.epsilon.Value{1};
isDebug     = isfield(sProcess.options, 'debug') && sProcess.options.debug.Value;
isAnatomy   = sProcess.options.rsspace.Value-1;
if isAnatomy==2
    isAnatomy = 0;
    isIdentity = 1;
else
    isIdentity = 0;
end
    
isApplySsp  = sProcess.options.applyssp.Value;

StudyName = sProcess.options.comment.Value;
if strcmp('@intra', StudyName)
    isNewStudy = 0;
else
    isNewStudy = 1;
end

% ===== ANALYZE DATABASE =====
ChannelFiles = {};
ChannelMats  = {};
CortexFiles  = {};
HeadModels  = {};
Kernels  = {};
MegInterp    = {};
MegInterpInv    = {};
isChanEqual  = [];
iInputSkip   = [];

% Check all the input files
for iInput = 1:length(sInputs)
    % No channel file: ignore
    if isempty(sInputs(iInput).ChannelFile)
        iInputSkip(end+1) = iInput;
        continue;
    end
    % Check channel file
    if ~any(file_compare(sInputs(iInput).ChannelFile, ChannelFiles))
        % Read channel file
        chanMat = in_bst_channel(sInputs(iInput).ChannelFile);
        % Check that same number of sensors
        nChan = length(chanMat.Channel);
        if ~isempty(ChannelMats) && (nChan ~= length(ChannelMats{1}.Channel))
            iInputSkip(end+1) = iInput;
        end
        % Check that for all the channel files have the same sizes for the Loc fields
        if ~isempty(ChannelMats) && ~isequal(cellfun(@(c)size(c),{chanMat.Channel.Loc},'UniformOutput',0), cellfun(@(c)size(c),{ChannelMats{1}.Channel.Loc},'UniformOutput',0))
            iInputSkip(end+1) = iInput;
        end
        % Check file types
        if ~any(ismember({'MEG','MEG GRAD','MEG MAG'}, unique({chanMat.Channel.Type})))
            %iInputSkip(end+1) = iInput;
            %continue;
        end
        % Remove head points
        if isfield(chanMat, 'HeadPoints')
            chanMat = rmfield(chanMat, 'HeadPoints');
        end
        % Add channel file to list
        ChannelFiles{end+1} = file_win2unix(sInputs(iInput).ChannelFile);
        ChannelMats{end+1}  = chanMat;
        % Check if channel description is the same as the first one
        if isempty(isChanEqual)
            isChanEqual = 1;
        else
            isChanEqual(end+1) = isequal({chanMat.Channel.Loc}, {ChannelMats{1}.Channel.Loc});
        end
        % Use anatomy as reference source space?
        if isAnatomy
            % find cortical surface file
            sSubj = load(file_fullpath(sInputs(iInput).SubjectFile), 'Cortex');
            CortexFiles{end+1} = sSubj.Cortex;
            sStudy = bst_get('Study', sInputs(iInput).iStudy);
            if length(sStudy.HeadModel)==1
                HeadModels{end+1} = sStudy.HeadModel.FileName;
            end
            tmp = sStudy.Result(~cat(1,sStudy.Result.isLink));
            for iR=1:length(tmp)
                if ~isempty(strfind(tmp(iR).FileName,'KERNEL'))
                    Kernels{end+1} = tmp(iR).FileName;
                    break
                end
            end
        else
            CortexFiles{end+1} = [];
        end
        % check bad channels
        sMat = in_bst_data(sInputs(1).FileName, 'ChannelFlag');
        channel_flag(:,iInput) = sMat.ChannelFlag;
    end
end

% Remove studies that cannot be processed
if ~isempty(iInputSkip)
    sInputs(iInputSkip) = [];
end
% Check that there is something to process
if isempty(sInputs)
    iFileChannel = []; AvgChannelMat = [];
    return;
elseif (length(ChannelFiles) == 1)
    iFileChannel = []; AvgChannelMat = [];
    return;
end

%% Average channel structure

% Base: first channel file in the list
AvgChannelMat = ChannelMats{1};
% Compute average channel structure (include ALL the sensor types)
if isAvgChan && ~all(isChanEqual)
    [AvgChannelMat, ~] = channel_average(ChannelMats);
    if isempty(AvgChannelMat)
%         iFileChannel = []; AvgChannelMat = [];
%         return;
        AvgChannelMat = ChannelMats{1};
    end
    % Consider that channel files  are all different from the average
    isChanEqual = 0 * isChanEqual;
end

%% ===== COMPUTE TRANSFORMATION =====

% Get MEG channels
iMeg = good_channel(AvgChannelMat.Channel, [], 'MEG');
% For each channel file: compute an interpolation matrix
for iFile = 1:length(ChannelMats)
    % If positions of sensors changed, or SSPs should be applied
    if ~isChanEqual(iFile) || isApplySsp || isDebug
        % Find the best-fitting sphere to source sensor locations
        chanlocs = figure_3d('GetChannelPositions', ChannelMats{iFile}, iMeg);
        [bfs_center, bfs_radius] = bst_bfs(chanlocs);
        % Compute combined SSP projector
        if isApplySsp && ~isempty(ChannelMats{iFile}.Projector)
            Proj = process_ssp2('BuildProjector',ChannelMats{iFile}.Projector,1);
            if ~isempty(Proj)
                Proj = Proj(iMeg,iMeg);
            else
                Proj = eye(numel(iMeg));
            end
        else
            Proj = [];
        end
        % Compute interpolation
        iMegSrc{iFile} = good_channel(ChannelMats{iFile}.Channel, [], 'MEG');
        [MegInterp{iFile}, MegInterpInv{iFile}, ~] = ComputeInterp(ChannelMats{iFile}.Channel(iMegSrc{iFile}), ...
            AvgChannelMat.Channel(iMeg), bfs_center, bfs_radius, [], epsilon, CortexFiles{iFile}, Proj);
    else
        MegInterp{iFile} = []; MegInterpInv{iFile} = [];
    end
end

%% ===== UPDATE DATA FILES =====
% Process each input data file
for iInput = 1:length(sInputs)
    % Get channel file that corresponds to this data file
    iFileChannel(iInput) = find(file_compare(sInputs(iInput).ChannelFile, ChannelFiles));
end

%% ===== SAVE CHANNEL FILE =====

% add coregistration operators to channel file
coreg = struct('isAvgChan', isAvgChan, ...
               'epsilon', epsilon, ...
               'isAnatomy', isAnatomy, ...
               'Study', struct('iStudy', [], 'MegInterp', [], ...
                               'MegInterpInv', [], 'ChannelFile', [], ...
                               'iMeg', []));
for iStudy = 1:numel(MegInterp)
    iFile = find(iFileChannel==iStudy,1);
    coreg.Study(iStudy).iStudy = sInputs(iFile).iStudy;
    if isIdentity
        coreg.Study(iStudy).MegInterp = eye(size(MegInterp{iStudy}));
        coreg.Study(iStudy).MegInterpInv = eye(size(MegInterpInv{iStudy}));
    else
        coreg.Study(iStudy).MegInterp = MegInterp{iStudy};
        coreg.Study(iStudy).MegInterpInv = MegInterpInv{iStudy};
    end
    coreg.Study(iStudy).ChannelFile = ChannelFiles{iStudy};
    coreg.Study(iStudy).iMeg = iMegSrc{iStudy};
end
AvgChannelMat.Coreg = coreg;

% save channel file
% make new condition if requested, otherwise save in @intra
if isApplySsp
    for iProj = 1:numel(AvgChannelMat.Projector)
        AvgChannelMat.Projector(iProj).Status = 0;
    end
    AvgChannelMat.Projector = [];
end
if isNewStudy
    iIntra = db_add_condition(sInputs(1).SubjectName, StudyName, 1);
else
    [~, iSubject] = bst_get('Subject', sInputs(1).SubjectName);
    [~, iIntra] = bst_get('AnalysisIntraStudy', iSubject);
end
% Add a history entry
if isAvgChan
    historyMsg = sprintf('New channels positions: average of %d channel files', length(chanMat));
else
    historyMsg = ['Using first channel file: ', sInputs(1).ChannelFile];
end
AvgChannelMat = bst_history('add', AvgChannelMat, 'meg_register', historyMsg);
db_set_channel(iIntra, AvgChannelMat, 2, 0);
OutputFiles = {sInputs.FileName};

end

function [WExtrap, WExtrapInv, src_xyz] = ComputeInterp(srcChan, destChan, bfs_center, bfs_radius, Whitener, epsilon, CortexFile, Proj)
    % ===== PARSE INPUTS =====
    if (nargin < 8) || isempty(Proj)
        Proj = []; 
    end    
    if (nargin < 7) || isempty(CortexFile)
        CortexFile = []; 
    end
    if (nargin < 6) || isempty(epsilon)
        epsilon = 0.0001;
    end
    if (nargin < 5) || isempty(Whitener)
        Whitener = [];
    end
    if (nargin < 4) || isempty(bfs_radius)
        bfs_radius = 0.07;
    end
    if (size(bfs_center,2) == 1)
        bfs_center = bfs_center';
    end
    WExtrap = [];

    % ===== COMPUTE SOURCE SPACE =====
    % If available, compute source grid from cortex anatomy, using
    % default options
    if ~isempty(CortexFile)
        Options.Method        = 'isotropic';
        Options.nLayers       = 5;
        Options.Reduction     = 3;
        Options.nVerticesInit = 500;
        Options.Resolution    = 0.01; % Isotropic option
        [sEnvelope, sCortex] = tess_envelope(CortexFile, 'convhull', Options.nVerticesInit, .001, []);
        src_xyz = bst_sourcegrid(Options, file_fullpath(CortexFile), sCortex, sEnvelope)';
    else
        % Computing spherical volume sourcespace.
        bfs_radius = 0.07;
        src_xyz = channel_extrapm('ComputeSphereGrid', bfs_center, bfs_radius)';
    end
    
    % ===== COMPUTE LEADFIELDS =====
    if ~isempty(CortexFile)
        % I go through the basic steps of forward modelling, maybe call to
        % bst_headmodeler instead?
        % Maybe build a single sphere model that spans the whole cortex,
        % but still within the sensors?
        % Create innerskull surface based on the cortex surface
        sSurfInner = tess_envelope(file_fullpath(CortexFile), 'convhull', 1082, .003);
        % Best fitting spheres (separate for both channel files..?)
        sph_orig = bst_os(srcChan, double(sSurfInner.Vertices), double(sSurfInner.Faces));
        sph_target = bst_os(destChan, double(sSurfInner.Vertices), double(sSurfInner.Faces));
        % Compute headmodels
        Gsrc2orig   = bst_meg_sph(src_xyz, srcChan,  sph_orig);
        Gsrc2target = bst_meg_sph(src_xyz, destChan, sph_target);
    else
        % Single sphere model
        Param = struct('Center', bfs_center', 'Radii', bfs_radius);
        % Compute headmodels
        Gsrc2orig   = bst_meg_sph(src_xyz, srcChan,  repmat(Param, [1,length(srcChan)]));
        Gsrc2target = bst_meg_sph(src_xyz, destChan, repmat(Param, [1,length(destChan)]));
    end
    if isempty(Gsrc2orig) || isempty(Gsrc2target)
        disp('EXTRAP> Error: One of the leadfields was not computed properly.');
        return
    end
    if ~isempty(Proj)
        Gsrc2orig = Proj * Gsrc2orig;
    end
    WExtrap = InterpFromHeadmodels(Gsrc2orig, Gsrc2target, Whitener, epsilon);
    WExtrapInv = InterpFromHeadmodels(Gsrc2target, Gsrc2orig, Whitener, epsilon);
end
    
function WExtrap = InterpFromHeadmodels(Gsrc2orig, Gsrc2target, Whitener, epsilon)
    % Apply whitener to leadfield
    if ~isempty(Whitener)
        Gsrc2orig = Whitener * Gsrc2orig;
    end

    % Computing SVD of Grammian and truncating based on epsilon.
    [U,S,V] = svd(Gsrc2orig * Gsrc2orig');
    s = diag(S);
    ss = sum(s);
    sss = 1 - (cumsum(s)./ss);
    i = find(sss<epsilon);
    i = i(1);
    si = 1./s;
    si = si(1:i);
    %display(['Using ' num2str(i) ' singular values/vectors (i.e, dimensions)']);

    % Computing extrapolation matrix.
    coefs = V(:,1:i) * diag(si) * U(:,1:i)';
    WExtrap = Gsrc2target * Gsrc2orig' * coefs;

    % Apply whitener to the interpolator
    if ~isempty(Whitener)
        WExtrap = WExtrap * Whitener;
    end
end


