function varargout = process_subscan_contributing_topographies( varargin )
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
% Authors: Peter Donhauser 2017

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
% Description the process
sProcess.Comment     = 'Subspace scan: Plot contributing topographies';
sProcess.FileTag     = '';
sProcess.Category    = 'File';
sProcess.SubGroup    = 'Sources';
sProcess.Index       = 113;
% Definition of the input accepted by this process
sProcess.InputTypes  = {'results'};
sProcess.OutputTypes = {''};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
sProcess.isSeparator = 0;
% Time point
sProcess.options.coord.Comment = 'MNI coordinates [x,y,z]:  ';
sProcess.options.coord.Type    = 'value';
sProcess.options.coord.Value   = {[0 0 0], 'mm', 1};     
end



%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFile = Run(sProcess, sInput)

coi = sProcess.options.coord.Value{1}./1000;

% load main files
sSubject = bst_get('Subject', sInput.SubjectFile);
sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);
sMat = in_bst_results(sInput.FileName);
nSubj = sMat.nAvg;

% find closest grid point
coi_scs = cs_convert(sMri, 'mni', 'scs', coi);
[~,iGrid] = min(sum(bsxfun(@minus, sMat.GridLoc, coi_scs).^2,2));

% load single subject subcorr maps
iHistory = find(strcmp(sMat.History(:,3), 'List of averaged files:'));
indHistory = iHistory+1:iHistory+sMat.nAvg;
for iSubj = 1:nSubj
    tmp = str_split(sMat.History{indHistory(iSubj),3}, ' ');
    sMatSubj = in_bst_results(tmp{2});
    sc(:,iSubj) = sMatSubj.ImageGridAmp(:,1);
    Options(iSubj) = sMatSubj.Options;
end

% Prepare subspace scanning
sProcessScan = process_scan_subspace('GetDescription');
sProcessScan.options.rsspace.Value = 1; % Set dimensions manually
sProcessScan.options.niter.Value{1} = 1;
sProcessScan.options.issave.Value = 0;

nRows = 5;
[~, indSubj] = sort(sc(iGrid,:));

close all

for iSubjSel = 1:numel(indSubj)
    iSubpl = rem(iSubjSel-1,nRows);
    if iSubpl==0
        fig = figure; clf; set(gcf,'Position', [680 0 560 950]);
    end
    iSubj = indSubj(iSubjSel);
    sProcessScan.options.ndim.Value{1} = Options(iSubj).nDim;
    if Options(iSubj).islast
        sProcessScan.options.order.Value = 2; % last
    else
        sProcessScan.options.order.Value = 1; % first
    end

    sInputScan = bst_process('GetInputStruct', Options(iSubj).SubspaceFile);
    SubjectNames{iSubjSel} = sInputScan.SubjectName;
    sSubscan = process_scan_subspace('Run', sProcessScan, sInputScan);
    [sc_coi,xx,yy] = subcorr(sSubscan.Options.G(:,:,iGrid), sSubscan.Options.Us, 1);
    g_coi = sSubscan.Options.Gorig(:,:,iGrid)*xx;
    
    hFig = view_topography(sInputScan.FileName,'MEG','2DDisc',g_coi);
    figure(fig)
    ha = subplot(nRows,2,1+iSubpl*2); axis off
    extract_topo_nice(hFig, ha);
    ylabel(sInputScan.SubjectName); 
    if iSubpl==0, title('Forward model'); end
    us_coi = sSubscan.Options.Usorig*yy;
    close(hFig)
    
    hFig = view_topography(sInputScan.FileName,'MEG','2DDisc',us_coi);
    figure(fig)
    ha = subplot(nRows,2,2+iSubpl*2); axis off
    extract_topo_nice(hFig, ha);
    if iSubpl==0, title('Subspace rotation'); end
    close(hFig)
    yl = ylabel(sprintf('Subcorr: %.2f',sc(iGrid,iSubj)));

end

OutputFile = sInput.FileName;
    
end

function ha = extract_topo_nice(hFig, ha)

posvals = get(hFig, 'Position');
set(hFig, 'Position', [posvals(1:2) 900 900]);
aa = get(hFig, 'Children');
set(aa,'XLim',[-1.5 1.5],'YLim',[-1.5 1.5]);
aaa = get(aa(end), 'Children');
c = 1;
for iChild = 1:numel(aaa)
    if strcmp(aaa(iChild).Type, 'line')
        set(aaa(iChild),'LineWidth', .01);
        topolines(c) = aaa(iChild);
        c = c + 1;
    end
end

frameGfx = getframe(hFig);
img = frameGfx.cdata;
axes(ha);
im = imshow(img);

hold on
factor1 = 342;
plot(-topolines(1).XData*factor1+430, -topolines(1).YData*factor1+450, 'k','LineWidth',2);
plot(-topolines(2).YData*factor1+430, -topolines(2).XData*factor1+450, 'k','LineWidth',2);
plot(-topolines(3).YData*factor1+430, -topolines(3).XData*factor1+450, 'k','LineWidth',2);
plot(-topolines(4).YData*factor1+430, -topolines(4).XData*factor1+450, 'k','LineWidth',2);

end












