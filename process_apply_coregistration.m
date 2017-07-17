function varargout = process_apply_coregistration( varargin )

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
    sProcess.Comment     = 'Apply co-registration';
    sProcess.FileTag     = '| coregistered';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'iES plugin';
    sProcess.Index       = 1001;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.processDim  = 1;   % Process channel by channel
    sProcess.Description = '';
    % Definition of the options
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
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

StudyName = sProcess.options.comment.Value;
if strcmp('@intra', StudyName)
    isNewStudy = 0;
    [~, iSubject] = bst_get('Subject', sInputs(1).SubjectName);
    [sIntra, iIntra] = bst_get('AnalysisIntraStudy', iSubject);    
else
    isNewStudy = 1;
    iIntra = db_add_condition(sInputs(1).SubjectName, StudyName, 1);
    sIntra = bst_get('Study', iIntra);
end

% Make sure data only come from one subject
Subjects = unique({sInputs.SubjectName});
if numel(Subjects)~=1
    bst_error('Multiple subjects not supported.');
end

% load @intra channel
iStudies = unique([sInputs.iStudy]);
% [~, iSubject] = bst_get('Subject', sInputs(1).SubjectName);
% [sIntra, iIntra] = bst_get('AnalysisIntraStudy', iSubject);
% ChannelFileIntra = sIntra.Channel.FileName;
ChannelFileIntra = sIntra.Channel.FileName;
ChannelMat = in_bst_channel(ChannelFileIntra);

% % make new condition if requested
% if isNewStudy
%     iIntra = db_add_condition(sInputs(1).SubjectName, StudyName, 1);
%     ChannelFileIntra = db_set_channel(iIntra, ChannelMat, 2, 0);
% end

% loop through the different studies
iOut = 1;
for iStudyCounter = 1:numel(iStudies)
    
    % find correct study indices
    iStudy = iStudies(iStudyCounter);
    sInputsStudy = sInputs([sInputs.iStudy]==iStudy);
    sInput = sInputsStudy(1);
    
    % check if co-registration was computed including this study
    if ~isfield(ChannelMat, 'Coreg') || ...
            ~(isempty(setdiff(iStudy,[ChannelMat.Coreg.Study.iStudy])))
        bst_error('Compute coregistration across all runs first');
    end
    
    % get correct co-registration operator
    [~,~,iStudyCoreg] = intersect(iStudy, [ChannelMat.Coreg.Study.iStudy]);
    MegInterp = ChannelMat.Coreg.Study(iStudyCoreg).MegInterp;
    iMeg = ChannelMat.Coreg.Study(iStudyCoreg).iMeg;
    
    % loop through all the files in this study
    for iFile = 1:numel(sInputsStudy)
        % load data and apply operator
        sMat = in_bst_data(sInputsStudy(iFile).FileName);
        sMat.F(iMeg,:) = MegInterp * sMat.F(iMeg,:);
        sMat.Comment = [sMat.Comment ' | coreg' num2str(sInput.iStudy)];
        sMat = bst_history('add', sMat, 'coregistered', [' - ' sInput.FileName]);
        
        % save data file in new study
        OutputFile = bst_process('GetNewFilename', bst_fileparts(ChannelFileIntra), 'data_coreg');
        % Output filename: add file tag
        OutputFile = file_unique(OutputFile);
        % Save file
        bst_save(OutputFile, sMat, 'v6');
        % Add file to database structure
        db_add_data(iIntra, OutputFile, sMat);
        OutputFiles{iOut} = OutputFile;
        iOut = iOut+1;
    end
    db_reload_studies(iIntra);
end

end