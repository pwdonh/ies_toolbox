function [X, iGoodChan, fs, indSegment] = get_data(Files, iGoodChan, TimeWindow, PassBand, WindowFcn, WinLength, WinOverlap)

% How many and which file type?
nFile = numel(Files);
[~, MatName] = in_bst(Files{1});
if strcmp(MatName,'F')
    ftype = 'data';
elseif strcmp(MatName,'Value')
    ftype = 'matrix';
end

% Get good channels (update on the fly, keep only the ones which are always good..)
for iFile = 1:nFile
    if strcmp(ftype,'data')
        sMat = in_bst_data(Files{iFile}, {'ChannelFlag', 'Time', 'Events'});
        iGoodChan = intersect(find(sMat.ChannelFlag == 1), iGoodChan)';
    else
        sMat = in_bst_matrix(Files{iFile});
    end
    % Get bad segments if any
    good{iFile} = ones(1, size(sMat.Time,2));
    if isfield(sMat,'Events')&&~isempty(sMat.Events)&&~isempty(sMat.Events(strcmp({sMat.Events.label}, 'BAD')))
        bad{iFile} = sMat.Events(strcmp({sMat.Events.label}, 'BAD')).samples;
        for iEp = 1:size(bad{iFile},2)
            if bad{iFile}(1,iEp)==0; bad{iFile}(1,iEp) = 1; end
            good{iFile}(bad{iFile}(1,iEp):bad{iFile}(2,iEp)) = 0;
        end
    else
        bad{iFile} = [];
    end
    if iFile==1
        iSample{1} = find((sMat.Time<TimeWindow(2))&(sMat.Time>TimeWindow(1)));
        t = sMat.Time;
        fs = round(1/(t(2)-t(1)));
    else
        if isempty(WinLength)&&(length(iSample{1}) ~= length(find((sMat.Time<TimeWindow(2))&(sMat.Time>TimeWindow(1)))))
            error('Make sure all trials have same time definition');
        end
        if fs ~= round(1/(sMat.Time(2)-sMat.Time(1)))
            error('Make sure all trials have same sampling rate');
        end
        iSample{iFile} = find((sMat.Time<TimeWindow(2))&(sMat.Time>TimeWindow(1)));
    end
end

% Window length and overlap
if ~isempty(WinLength)
    Lwin  = round(WinLength * fs);
    Loverlap = round(Lwin * WinOverlap / 100);
else
    Lwin = numel(iSample{1});
end

% Compute window function
if ~isempty(WindowFcn)
    win = bst_window(WindowFcn, Lwin)';
else
    win = ones(1,Lwin);
end

% Load data
X = zeros(numel(iGoodChan), Lwin, nFile); %%% todo: compute number of segments here
nBad = 0;
iSegment = 1;
for iFile = 1:numel(Files)
    [sMat, MatName] = in_bst(Files{iFile});
    Xtmp = sMat.(MatName)(iGoodChan, iSample{iFile});
    if ~isempty(WinLength)
        nTime = size(Xtmp,2);
        nWin = floor((nTime - Loverlap) ./ (Lwin - Loverlap));
        for iWin = 1:nWin
            iTimes = (1:Lwin) + (iWin-1)*(Lwin - Loverlap);
            % Check if this segment is outside of ALL the bad segments (either entirely before or entirely after)
            if ~isempty(bad{iFile}) && (~all((iTimes(end) < bad{iFile}(1,:)) | (iTimes(1) > bad{iFile}(2,:))))
                disp(sprintf('BST> Skipping window #%d because it contains a bad segment.', iWin));
                nBad = nBad + 1;
                continue;
            end
            % subtract mean before windowing
            Xtmp2 = bsxfun(@minus, Xtmp(:,iTimes), mean(Xtmp(:,iTimes),2));
            X(:,:,iSegment) = bst_bsxfun(@times, Xtmp2, win);
            iSegment = iSegment + 1;
        end
        % make index for files and windows
        tmpindex = iSegment-nWin+nBad:iSegment-1;
        indSegment(tmpindex,:) = [1:nWin-nBad; repmat(iFile,1,nWin-nBad)]';
    else
        X(:,:,iFile) = bst_bsxfun(@times, Xtmp, win);
        indSegment(iFile,:) = [1 iFile];
    end
end

if ~isempty(PassBand)
    for iSegment = 1:size(X,3)
        X(:,:,iSegment) = process_bandpass('Compute', X(:,:,iSegment), fs, PassBand(1), PassBand(2), 'bst-fft-fir', 1);
    end
end

% Zero-out bad segments..
if isempty(WinLength)
    for iFile = 1:numel(Files)
        X(:,good{iFile}(iSample{iFile})<1,iFile) = 0;
    end
end

end