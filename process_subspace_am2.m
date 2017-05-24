function varargout = process_subspace_coh2( varargin )
%
%
%

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
% Authors: Peter Donhauser, 2015

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Compute Signal Subspace: Amplitude Modulation';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Connectivity';
    sProcess.Index       = 651;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'matrix','data'};
    sProcess.OutputTypes = {'data'};
    sProcess.nInputs     = 2;
    sProcess.nMinFiles   = 1;
    sProcess.isPaired    = 1;
    
    % === TITLE
    sProcess.options.label2.Comment = '<BR><U><B>Estimator options</B></U>:';
    sProcess.options.label2.Type    = 'label';
    % Time window
    sProcess.options.timewindow.Comment = 'Time window of interest: ';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];    
    % Row index
    sProcess.options.irow.Comment = 'Row index:  ';
    sProcess.options.irow.Type    = 'value';
    sProcess.options.irow.Value   = {1, '', 0}; % what do the cell entries have to be?       
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
    sProcess.options.isoutf.Comment = 'Output filter files';
    sProcess.options.isoutf.Type    = 'checkbox';
    sProcess.options.isoutf.Value   = 0;    
    % 
    sProcess.options.issavecov.Comment = 'Save single-trial covariance matrices';
    sProcess.options.issavecov.Type    = 'checkbox';
    sProcess.options.issavecov.Value   = 0;        
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputA, sInputB) %#ok<DEFNU>

    % Get inputs
    sFilesRef = {sInputA.FileName}; % reference signals (from a matrix file)
    sFilesDat = {sInputB.FileName}; % data (MEG signals)
    if numel(sFilesDat)~=numel(sFilesDat)
        bst_error('Needs same number of inputs.');
    end

    Options.Stimulus  = sProcess.options.timewindow.Value{1};
    Options.PassBand  = sProcess.options.bandpass.Value{1};   
    Options.alpha     = sProcess.options.reg.Value{1};
    Options.epsilon   = sProcess.options.epsilon.Value{1};
    Options.isOutF    = sProcess.options.isoutf.Value;
    Options.WindowFcn = sProcess.options.wfcn.Value;
    Options.isSaveCov = sProcess.options.issavecov.Value;
	isSplit = sProcess.options.split.Value;
    if isSplit
        Options.WinLength = sProcess.options.win_length.Value{1};
        Options.WinOverlap = sProcess.options.win_overlap.Value{1};
    else
        Options.WinLength = [];
        Options.WinOverlap = [];
    end    
    Options.iRow      = sProcess.options.irow.Value{1};
    Options.isAM      = 1;
    Options.isSeedAM  = 0;
    Options.isProj    = 0;
    
    % Call function
    OutputFiles = compute_signal_subspace(sFilesDat, sFilesRef, Options);
    
end



