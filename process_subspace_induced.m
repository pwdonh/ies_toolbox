function varargout = process_subspace_induced( varargin )
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
    sProcess.Comment     = 'Compute Signal Subspace: Induced responses';
    sProcess.FileTag     = '| subsp';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'iES plugin';
    sProcess.Index       = 100;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
%     sProcess.isSeparator = 1;
    % Filter
    sProcess.options.bandpass.Comment = 'Frequency band: ';
    sProcess.options.bandpass.Type    = 'range';
    sProcess.options.bandpass.Value   = {[8, 13], 'Hz', 2};
    % Stimulus time window
    sProcess.options.stim.Comment = 'Time window of interest: ';
    sProcess.options.stim.Type    = 'range';
    sProcess.options.stim.Value   = {[0, 2], 's', []};
    % Baseline time window
    sProcess.options.base.Comment = 'Reference time window:';
    sProcess.options.base.Type    = 'range';
    sProcess.options.base.Value   = {[-2, 0], 's', []};
    % PCA
    sProcess.options.ispca.Comment = 'Use PCA instead (no reference)';
    sProcess.options.ispca.Type    = 'checkbox';
    sProcess.options.ispca.Value   = 0; 
    % Regularization parameter
    sProcess.options.reg.Comment = 'Covariance regularization:  ';
    sProcess.options.reg.Type    = 'value';
    sProcess.options.reg.Value   = {.05, '', 3};
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
    sProcess.options.wfcn.Value   = 'tukey';
    % Save spatial filters (by default only spatial patterns)
    sProcess.options.isoutf.Comment = 'Save spatial filter files';
    sProcess.options.isoutf.Type    = 'checkbox';
    sProcess.options.isoutf.Value   = 0;    
    % Save all covariance matrices; those are needed for permutation tests
    sProcess.options.issavecov.Comment = 'Save single-trial covariance matrices (for statistics)';
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
    isPca = sProcess.options.ispca.Value;
    if isPca
        Options.Baseline = [];
    else
        Options.Baseline  = sProcess.options.base.Value{1};
    end
	isSplit = sProcess.options.split.Value;
    if isSplit
        Options.WinLength = sProcess.options.win_length.Value{1};
        Options.WinOverlap = sProcess.options.win_overlap.Value{1};
    else
        Options.WinLength = [];
        Options.WinOverlap = [];
    end
    Options.Stimulus  = sProcess.options.stim.Value{1};
    Options.alpha     = sProcess.options.reg.Value{1};
    Options.WindowFcn = sProcess.options.wfcn.Value;
    Options.isOutF    = sProcess.options.isoutf.Value;
    Options.isSaveCov = sProcess.options.issavecov.Value;
    Options.isProj    = 0;
    
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
        OutputFiles = compute_signal_subspace(sFiles, [], Options);
    end
    
end


