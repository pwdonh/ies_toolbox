function varargout = process_apply_spfilter2( varargin )
% PROCESS_SPATIAL_FILTER2: Apply one or several spatial filters to a sensor
% data file. File A contains filters, File B the data.
%
% USAGE:  OutputFiles = process_compare2('Run', sProcess, sInputsA, sInputsB)

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
% Authors: Peter Donhauser, 2015

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
    % Description the process
    sProcess.Comment     = 'Apply spatial filter';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'iES plugin';
    sProcess.Index       = 302;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'matrix'};
    sProcess.nInputs     = 2;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    
    % === Number of components
    sProcess.options.ncomp.Comment = 'Number of components:  ';
    sProcess.options.ncomp.Type    = 'value';
    sProcess.options.ncomp.Value   = {1, '', 0};
    % === Order of components
    sProcess.options.label4.Type    = 'label';
    sProcess.options.label4.Comment = '<HTML><BR>Order :';
    sProcess.options.order.Comment = {'first', 'last'};
    sProcess.options.order.Type    = 'radio';
    sProcess.options.order.Value   = 0;    
    % === Matrix ordering
    sProcess.options.matrix.Type    = 'label';
    sProcess.options.matrix.Comment = '<HTML><BR>Rows in output represent :';
    sProcess.options.matrix.Comment = {'Epochs (one matrix per component)', ...
                                        'Components (one matrix per epoch)'};
    sProcess.options.matrix.Type    = 'radio';
    sProcess.options.matrix.Value   = 0;
    % === Comment
    sProcess.options.comment.Comment = 'Comment:';
    sProcess.options.comment.Type    = 'text';
    sProcess.options.comment.Value   = '';    
    % Normalize
    sProcess.options.isnoisenorm.Comment = 'Normalize values using noise statistics';
    sProcess.options.isnoisenorm.Type    = 'checkbox';
    sProcess.options.isnoisenorm.Value   = 0;    
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputsA, sInputsB)

Options.nComp = sProcess.options.ncomp.Value{1};
Options.Comment = sProcess.options.comment.Value;
Options.order = sProcess.options.order.Value;
Options.matrix = sProcess.options.matrix.Value;
isNoiseNorm = sProcess.options.isnoisenorm.Value;

if isempty(Options.Comment)
    Options.Comment = 'Spatial filter';
end

if isNoiseNorm
    % The noise covariance matrix is used to normalize the signal power
    sStudy = bst_get('Study', sInputsA.iStudy);
    NoiseFile = sStudy.NoiseCov.FileName;
else
    NoiseFile = [];
end

OutputFiles = apply_spatial_filter(sInputsA.FileName, ...
    {sInputsB.FileName}, NoiseFile, Options);
% I need to output this, to avoid error in bst_process()
% OutputFiles = {sInputsB.FileName};
% db_reload_conditions(sInputsB(1).iStudy);

end

function MatrixFiles = apply_spatial_filter(FilterFile, DataFiles, NoiseFile, Options)

sInput = bst_process('GetInputStruct',DataFiles{1});
sInputA = bst_process('GetInputStruct',FilterFile);
% Load filter file
sMat = in_bst_data(FilterFile);
% Check that it's an appropriate file
if ~isfield(sMat, 'C1')||~isfield(sMat, 'GoodChannel')
    bst_error('Not an appropriate subspace or dipole file');
end
if sInputA.iStudy~=sInput(1).iStudy
    ChannelMatFilt = in_bst_channel(sInputA.ChannelFile);
    ChannelMatData = in_bst_channel(sInput(1).ChannelFile);
    ChannelNames = {ChannelMatFilt.Channel(sMat.GoodChannel).Name};
    iChanDat = channel_find(ChannelMatData.Channel, ChannelNames);
    if length(sMat.GoodChannel)~=length(iChanDat)
        ChannelNames = {ChannelMatData.Channel(iChanDat).Name};
        iChanFilt = channel_find(ChannelMatFilt.Channel, ChannelNames);
        [~,~,iChanChan]=intersect(iChanFilt,sMat.GoodChannel);
    else
        iChanFilt = sMat.GoodChannel;
        iChanChan = 1:length(iChanFilt);
    end
else
    iChanDat = sMat.GoodChannel;
    iChanFilt = sMat.GoodChannel;
    iChanChan = 1:length(iChanFilt);
end

%% Compute filter kernel

% Get relevant topographies and noise covariance
if Options.order==1
    P = sMat.F(iChanFilt, 1:Options.nComp);
    Cn = sMat.C2(iChanChan,iChanChan);
elseif Options.order==2
    P = sMat.F(iChanFilt, end-Options.nComp+1:end);
    Cn = sMat.C1(iChanChan,iChanChan);
end
% Compute whitener
[U, S] = svd(Cn,'econ');
% regularize
if isfield(sMat.Options, 'epsilon')
    epsilon = sMat.Options.epsilon;
else
    epsilon = .05;
end
s = diag(S);
ss = sum(s);
sss = 1 - (cumsum(s)./ss);
i = find(sss<epsilon);
if isempty(i)
    i = length(sss);
else
    i = i(1);
end
Proj = sqrt(inv(S(1:i,1:i)))*U(:,1:i)';
% Compute filter
W = (Proj'*pinv(Proj*P)')';

%% Load noise covariance

if ~isempty(NoiseFile)
    NoiseMat = load(file_fullpath(NoiseFile));
    Nm = NoiseMat.NoiseCov(iChanFilt,iChanFilt);
    noisenorm = sqrt(diag(W*Nm*W'));
end

%% Apply to data files

for iFile = 1:numel(DataFiles)
    DataMat = in_bst_data(DataFiles{iFile});
    Events{iFile} = DataMat.Events;
    X = DataMat.F(iChanDat,:);
    if ~isempty(NoiseFile)
        if Options.matrix==1
            Y(:,:,iFile) = bsxfun(@rdivide, W * X, noisenorm);
        else
            Y{iFile} = bsxfun(@rdivide, W * X, noisenorm);
            t{iFile} = DataMat.Time;
        end
    else
        Y{iFile} = W * X;
        if Options.matrix==1
            Y(:,:,iFile) = W * X;
        else
            Y{iFile} = W * X;
            t{iFile} = DataMat.Time;
        end        
    end
end

%% Save resulting matrix files

if Options.matrix==1
    for iComp = 1:Options.nComp
        MatrixMat = db_template('matrixmat');
        MatrixMat.Value = squeeze(Y(iComp,:,:))';
        MatrixMat.Description = {};
        %     MatrixMat.Events = Events;
        for iFile = 1:numel(DataFiles)
            MatrixMat.Description{iFile,1} = ['Epoch ' num2str(iFile)];
        end
        MatrixMat.Time = DataMat.Time;
        MatrixMat.Comment = [Options.Comment ' component #' num2str(iComp)];
        % Get filename
        MatrixFiles{iComp} = bst_process('GetNewFilename', bst_fileparts(DataFiles{1}), 'matrix');
        % Output filename: add file tag
        MatrixFiles{iComp} = file_unique(MatrixFiles{iComp});
        % Save file
        bst_save(MatrixFiles{iComp}, MatrixMat, 'v6');
        db_add_data(sInput.iStudy, MatrixFiles{iComp}, MatrixMat);        
    end
elseif Options.matrix==2
    MatrixMat = db_template('matrixmat');
    MatrixMat.Description = {};
    for iComp = 1:Options.nComp
        MatrixMat.Description{iComp,1} = ['Component ' num2str(iComp)];
    end
    for iFile = 1:numel(DataFiles)
        MatrixMat.Value = Y{iFile};
        MatrixMat.Time = t{iFile};
        MatrixMat.Comment = [Options.Comment ' epoch (#' num2str(iFile) ')'];
        MatrixMat.Events = Events{iFile};
        % Get filename
        MatrixFiles{iFile} = bst_process('GetNewFilename', bst_fileparts(DataFiles{1}), 'matrix');
        % Output filename: add file tag
        MatrixFiles{iFile} = file_unique(MatrixFiles{iFile});
        % Save file
        bst_save(MatrixFiles{iFile}, MatrixMat, 'v6');
        db_add_data(sInput.iStudy, MatrixFiles{iFile}, MatrixMat);
    end
end

end