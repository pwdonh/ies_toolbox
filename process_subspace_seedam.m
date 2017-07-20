function varargout = process_subspace_stim( varargin )
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
    sProcess.Comment     = 'Compute Signal Subspace: Seed AM';
    sProcess.FileTag     = '| subsp';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'iES plugin';
    sProcess.Index       = 112;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Time window
    sProcess.options.timewindow.Comment = 'Time window of interest: ';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
    % Mni coordinate
    sProcess.options.coord.Comment = 'MNI coordinates [x,y,z]:  ';
    sProcess.options.coord.Type    = 'value';
    sProcess.options.coord.Value   = {[-39, -27, 55], 'mm', 1};
    % Use Dipole orientation
    sProcess.options.isestimate.Comment = 'Estimate dipole orientation from data';
    sProcess.options.isestimate.Type    = 'checkbox';
    sProcess.options.isestimate.Value   = 1;       
    % Set dipole orientation
    sProcess.options.orient.Comment = 'If not: dipole orientation [x,y,z]:  ';
    sProcess.options.orient.Type    = 'value';
    sProcess.options.orient.Value   = {[0 0 0], '', 1};
    % Filter
    sProcess.options.bandpass.Comment = 'Frequency band: ';
    sProcess.options.bandpass.Type    = 'range';
    sProcess.options.bandpass.Value   = {[13, 40], 'Hz', 2};
    % Regularization parameter
    sProcess.options.reg.Comment = 'Covariance regularization:  ';
    sProcess.options.reg.Type    = 'value';
    sProcess.options.reg.Value   = {0, '', 3}; % what do the cell entries have to be?
    % Regularization parameter
    sProcess.options.epsilon.Comment = 'Dimensionality reduction (epsilon):  ';
    sProcess.options.epsilon.Type    = 'value';
    sProcess.options.epsilon.Value   = {0.001, '', 3}; % what do the cell entries have to be?    
    % Split up
    sProcess.options.split.Comment = 'Split up into time windows';
    sProcess.options.split.Type    = 'checkbox';
    sProcess.options.split.Value   = 0;
    % Option: Window (Length)
    sProcess.options.win_length.Comment    = 'Window length: ';
    sProcess.options.win_length.Type       = 'value';
    sProcess.options.win_length.Value      = {1, 's', []};
    % Option: Window (Overlapping ratio)
    sProcess.options.win_overlap.Comment    = 'Window overlap ratio: ';
    sProcess.options.win_overlap.Type       = 'value';
    sProcess.options.win_overlap.Value      = {50, '%', 1};    
    % Window function
    sProcess.options.wfcn.Comment = 'Window function (hann, tukey,..):';
    sProcess.options.wfcn.Type    = 'text';
    sProcess.options.wfcn.Value   = 'hann';
    % 
    sProcess.options.isproj.Comment = 'Project out seed topography';
    sProcess.options.isproj.Type    = 'checkbox';
    sProcess.options.isproj.Value   = 0;        
    % 
    sProcess.options.isoutf.Comment = 'Output filter files';
    sProcess.options.isoutf.Type    = 'checkbox';
    sProcess.options.isoutf.Value   = 0;    
    % 
    sProcess.options.issavecov.Comment = 'Save single-trial covariance matrices';
    sProcess.options.issavecov.Type    = 'checkbox';
    sProcess.options.issavecov.Value   = 0;    
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs)

    % Get inputs
    sFiles = {sInputs.FileName};
    Options.PassBand  = sProcess.options.bandpass.Value{1};
    Options.NoiseBand  = [];
    Options.Stimulus  = sProcess.options.timewindow.Value;
    Options.Baseline = [];
	isSplit = sProcess.options.split.Value;
    if isSplit
        Options.WinLength = sProcess.options.win_length.Value{1};
        Options.WinOverlap = sProcess.options.win_overlap.Value{1};
    else
        Options.WinLength = [];
        Options.WinOverlap = [];
    end
    Options.alpha     = sProcess.options.reg.Value{1};
    Options.epsilon   = sProcess.options.epsilon.Value{1};
    Options.WindowFcn = sProcess.options.wfcn.Value;
    Options.isProj    = sProcess.options.isproj.Value;
    Options.isOutF    = sProcess.options.isoutf.Value;
    Options.isSaveCov = sProcess.options.issavecov.Value;
    Options.isAM      = 1;
    Options.isSeedAM  = 1;
    
    % Extract scout time-series and save in matrix files
    if ~sProcess.options.isestimate.Value
        orient = sProcess.options.orient.Value{1};
    else
        orient = [];
    end
    [sFilesRef, topo, ~, iVertex] = ExtractScoutTs(sProcess.options.coord.Value{1}, ...
        sInputs, Options.PassBand, Options.Stimulus, 1, orient);

    Options.proj = topo;
    Options.iRow = 1;
    Options.iVertex = iVertex;
    
    % Get data file names
    mint = 0; maxt = 1;
    for iInput = 1:numel(sInputs)
        tmp = str_split(sInputs(iInput).FileName,'|');
        sFilesDat{iInput} = tmp{3};
        t = in_bst(sFilesDat{iInput},'Time');
        mint = min([t mint]);
        maxt = max([t maxt]);
    end
    if isempty(Options.Stimulus)
        Options.Stimulus = [mint maxt];
    end
    
    % Handle multiple subjects
    [SubjectNames,~,subject_index] = unique({sInputs.SubjectName});
    if numel(SubjectNames)>1
        fprintf('Multiple subjects. Processing one by one..\n');
        OutputFiles = {};
        for iSubj = 1:numel(SubjectNames)
            fprintf('Subject %s\n', SubjectNames{iSubj});
            sFilesSubj = sFiles(subject_index==iSubj);
            % Call function
            OutputFilesSubj{iSubj} = compute_signal_subspace(sFilesSubj, [], Options);
        end
        OutputFiles = cat(2,OutputFilesSubj{:});
    else
        % Call function
        OutputFiles = compute_signal_subspace(sFilesDat, sFilesRef, Options);
    end
    
    % Delete reference files
    sInputsRef = bst_process('GetInputStruct', sFilesRef);
    sProcessDelete = process_delete('GetDescription');
    sProcess.options = struct_copy_fields(sProcess.options, sProcessDelete.options, 0);
    tmp = process_delete('Run', sProcess, sInputsRef);
        
end

function [output, topo, t, iVertex] = ExtractScoutTs(seedInfo, sInputs, BandPass, TimeWindow, isSave, orient)

if iscell(seedInfo)
    isMni = 0;
    AtlasName = sProcess.options.scouts.Value{1};
    ScoutName = sProcess.options.scouts.Value{2};
    if numel(ScoutName)>1
        error('No!');
    end    
else
    if length(seedInfo)==3
        isMni = 1;
    end
end

ResultMat = in_bst_results(sInputs(1).FileName);
SurfaceMat = in_tess_bst(ResultMat.SurfaceFile);
HeadModelMat = in_bst_headmodel(ResultMat.HeadModelFile);
ChannelMat = in_bst_channel(sInputs(1).ChannelFile);

if ~isMni
    iAtlas = strcmp({SurfaceMat.Atlas.Name}, AtlasName);
    iScout = strcmp({SurfaceMat.Atlas(iAtlas).Scouts.Label}, ScoutName{1});
    iVertex = SurfaceMat.Atlas(iAtlas).Scouts(iScout).Vertices;
    label = SurfaceMat.Atlas(iAtlas).Scouts(iScout).Label;
else
    sSubject = bst_get('Subject',sInputs.SubjectFile);
    sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);
    % find closest grid point
    coi_scs = cs_convert(sMri, 'mni', 'scs', seedInfo./1000);
    [~,iVertex] = min(sum(bsxfun(@minus, ResultMat.GridLoc, coi_scs).^2,2));
    label = sprintf('MNI [%d, %d, %d]', seedInfo);
end
rowindex = iVertex*3-[2 1 0];
iChan = ResultMat.GoodChannel;
% iChan = 1:numel(ChannelMat.Channel);

% forward and inverse fields
G = HeadModelMat.Gain(:,rowindex);
K = ResultMat.ImagingKernel(rowindex,:);

fprintf('Loading seed time-series, trial 000\n');
yf_all = [];
for iInput = 1:numel(sInputs)
    fprintf('\b\b\b\b%3d\n',iInput)
    tmp = str_split(sInputs(iInput).FileName,'|');
    DataMat = in_bst_data(tmp{3});
    comments{iInput} = DataMat.Comment;
    events{iInput} = DataMat.Events;
    t = DataMat.Time;
    timevec{iInput} = t;
    fs = round(1/(t(2)-t(1)));
    X = DataMat.F(iChan,:);
    Xf = process_bandpass('Compute', X, fs, BandPass(1), BandPass(2));
    if ~isempty(TimeWindow)
        tindex = (t<TimeWindow(2))&(t>TimeWindow(1));
    else
        tindex = 1:size(X,2);
    end
    y{iInput} = K*X;
    yf_all = [yf_all K*Xf(:,tindex)];
end
% yf_all = reshape(yf,3,size(yf,2)*numel(sInputs));
% y_all = reshape(y,3,size(y,2)*numel(sInputs));
[U,S,V] = svd(yf_all, 'econ');

% See if dipole orientation is given as input
if ~isempty(orient)
    U(:,1) = orient;
end

topo = G*U(:,1);
for iInput = 1:numel(sInputs)
    output{iInput} = U(:,1)'*y{iInput};
end

if isSave
    iStudy = unique([sInputs.iStudy]);
    MatrixMat = db_template('matrixmat');
    if ~isMni
        MatrixMat.Description = {[AtlasName ' - ' ScoutName{1}]};
    else
        MatrixMat.Description = {num2str(seedInfo)};
    end
    for iInput = 1:numel(sInputs)
        MatrixMat.Time = timevec{iInput};
        MatrixMat.Value = output{iInput};
        MatrixMat.Comment = [comments{iInput} ' | ' label];
        MatrixMat.Events = events{iInput};
        % Get filename
        MatrixFiles{iInput} = bst_process('GetNewFilename', bst_fileparts(tmp{3}), 'matrix');
        % Output filename: add file tag
        MatrixFiles{iInput} = file_unique(MatrixFiles{iInput});
        % Save file
        bst_save(MatrixFiles{iInput}, MatrixMat, 'v6');
        db_add_data(sInputs(iInput).iStudy, MatrixFiles{iInput}, MatrixMat);
    end
    output = file_short(MatrixFiles);
end    

end













