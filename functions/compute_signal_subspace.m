function [output, X1, X2, y] = compute_signal_subspace(FilesA, FilesB, Options)
% Computes the signal subspace for use in scanning procedures like MUSIC.
% Includes the standard SVD/PCA approach, but also contrasts between
% stimulus/baseline, different conditions and coherence with a reference
% signal (all solved using generalized eigenvalue decomposition).
%
% Peter Donhauser

%% === PREPARE

if ~isfield(Options, 'isSave')
    Options.isSave = 1;    
end

% Check which analysis type to use
if isempty(FilesB)
    Files = FilesA;
    if ~isempty(Options.Baseline)&&isempty(Options.NoiseBand)
        pipeline = 'stim_base';
    elseif ~isempty(Options.NoiseBand)&&isempty(Options.Baseline)
        pipeline = 'freq1_freq2';
    elseif isempty(Options.Baseline)
        pipeline = 'pca';        
    end
else
    [~, MatName] = in_bst(FilesB{1});
    if strcmp(MatName, 'F')
        Files = cat(2,FilesA,FilesB);
        file_index = [ones(1,numel(FilesA)) ones(1,numel(FilesB))*2];
        pipeline = 'cond1_cond2';
    elseif strcmp(MatName, 'Value')
        Files = FilesA; % data is processed similarly
        if ~Options.isAM
            pipeline = 'coh_ref';
        else
            if Options.isSeedAM
                pipeline = 'aec';
            else
                pipeline = 'coh_ref_am';
            end
        end
    end
end

% Make sure files only come from one (possibly coregistered) condition
sInputs = bst_process('GetInputStruct', Files);
iStudy = unique([sInputs.iStudy]);
if numel(iStudy)~=1
    bst_error('Input must be from one condition');
end
sStudy = bst_get('Study', sInputs.iStudy);

% Load additional files
ChannelMat = in_bst_channel(sInputs(1).ChannelFile);

% Get indices to relevant channels
iChan = good_channel(ChannelMat.Channel, [], 'MEG');
iGoodChan = iChan;

% Load data (+ apply windowing and filtering)
fprintf('Load data and filter the requested time windows\n');
switch pipeline
    case 'stim_base'
        [X1, iGoodChan, fs, Options.indSegment] = get_data(Files, iGoodChan, Options.Stimulus, Options.PassBand, Options.WindowFcn, Options.WinLength, Options.WinOverlap);
        X2 = get_data(Files, iGoodChan, Options.Baseline, Options.PassBand, Options.WindowFcn, Options.WinLength, Options.WinOverlap);
        % Cut to same size
        if size(X1,2)>size(X2,2)
            warning('Cutting the two time windows to same size');
            X1 = X1(:,1:size(X2,2),:);
        elseif size(X1,2)>size(X2,2)
            warning('Cutting the two time windows to same size');
            X2 = X2(:,1:size(X1,2),:);
        end
        y = [];
    case 'cond1_cond2'
        [X, iGoodChan, fs, Options.indSegment] = get_data(Files, iGoodChan, Options.Stimulus, Options.PassBand, Options.WindowFcn, Options.WinLength, Options.WinOverlap);
        X1 = X(:,:,file_index==1);
        X2 = X(:,:,file_index==2);
        y = [];
    case 'freq1_freq2'
        [X1, iGoodChan, fs, Options.indSegment] = get_data(Files, iGoodChan, Options.Stimulus, Options.PassBand, Options.WindowFcn, Options.WinLength, Options.WinOverlap);
        X2 = get_data(Files, iGoodChan, Options.Stimulus, Options.NoiseBand, Options.WindowFcn, Options.WinLength, Options.WinOverlap);
        y = [];
    case 'coh_ref'
        % Load data (+ apply windowing and filtering)
        [X1, iGoodChan, fs, Options.indSegment] = get_data(Files, iGoodChan, Options.Stimulus, [], Options.WindowFcn, Options.WinLength, Options.WinOverlap);
        % Load reference signals
        X2 = [];
        y = get_data(FilesB, Options.iRow, Options.Stimulus, [], Options.WindowFcn, Options.WinLength, Options.WinOverlap);
    case 'aec'
        % Load data (+ apply windowing and filtering)
        [X1, iGoodChan, fs, Options.indSegment] = get_data(Files, iGoodChan, Options.Stimulus, Options.PassBand, Options.WindowFcn, Options.WinLength, Options.WinOverlap);
        % Load reference signals
        X2 = [];
        y = get_data(FilesB, Options.iRow, Options.Stimulus, Options.PassBand, Options.WindowFcn, Options.WinLength, Options.WinOverlap);
        for iSegment = 1:size(y,3)
            y(1,:,iSegment) = abs(hilbert(y(1,:,iSegment)')');
        end
    case 'coh_ref_am'
        % Load data (+ apply windowing and filtering)
        [X1, iGoodChan, fs, Options.indSegment] = get_data(Files, iGoodChan, Options.Stimulus, Options.PassBandAM, Options.WindowFcn, Options.WinLength, Options.WinOverlap);
        X2 = [];
        % Load reference signals
        y = get_data(FilesB, Options.iRow, Options.Stimulus, Options.PassBand, Options.WindowFcn, Options.WinLength, Options.WinOverlap);
    case 'pca'
        % Load data (+ apply windowing and filtering)
        [X1, iGoodChan, fs, Options.indSegment] = get_data(Files, iGoodChan, Options.Stimulus, Options.PassBand, Options.WindowFcn, Options.WinLength, Options.WinOverlap);
        y = [];
        if ~isempty(sStudy.NoiseCov)
            NoiseCovMat = load(file_fullpath(sStudy.NoiseCov.FileName));
            Options.NoiseCov = NoiseCovMat.NoiseCov(iGoodChan,iGoodChan);
        else
            Options.NoiseCov = eye(size(X,1));
        end
end
Options.fs = fs;
if isfield(Options, 'proj')&&~isempty(Options.proj)&&...
        length(Options.proj)==numel(ChannelMat.Channel)
    Options.proj = Options.proj(iGoodChan,1);
end
fprintf('Ignoring %d bad channels.\n', numel(iChan)-numel(iGoodChan));

%% === COMPUTE

[SNR, Results, C1, C2, C1_all, C2_all, Pow_all, Pow_all_unreg, C_orig, y_all] = ...
    compute_signal_subspace_core(X1, X2, y, pipeline, Options);

%% === SAVE RESULTS

Comment = {'Subspace patterns', 'Subspace filters'};
FilesOut = {};
for iOutput = 1:Options.isOutF+1 % for filters and patterns
    % Create data structure
    DataMat = in_bst_data(Files{1});
    DataMat.F               = zeros(size(DataMat.F,1), size(Results{iOutput}, 2));
    DataMat.F(iGoodChan,:)  = Results{iOutput};
    DataMat.Comment         = [Comment{iOutput} ' - ' pipeline];
    DataMat.Time            = 1:size(Results{iOutput}, 2);
    DataMat.DataType        = 'recordings';
    DataMat.Device          = 'subspace';
    DataMat.SNR             = SNR;
    DataMat.C1              = C1;
    DataMat.C2              = C2;
    if Options.isSaveCov
        if strcmp(pipeline, 'coh_ref_am')||strcmp(pipeline, 'aec')
            DataMat.C_orig  = C_orig;
            DataMat.y_all   = y_all;
        else
            DataMat.C1_all      = C1_all;
            DataMat.C2_all      = C2_all;
        end
    end
    DataMat.Pow_all     = Pow_all;
    DataMat.Pow_all_unreg     = Pow_all_unreg;
    DataMat.GoodChannel     = iGoodChan; % put this field like in an imaging kernel
    DataMat.Pipeline        = pipeline;
    DataMat.Options         = Options;
    DataMat.FilesA          = FilesA;
    DataMat.FilesB          = FilesB;
    DataMat.ChannelFlag     = ones(size(DataMat.ChannelFlag))*-1;
    DataMat.ChannelFlag(iGoodChan) = 1;
    % Save in database
    if Options.isSave
        OutputFile = bst_process('GetNewFilename', bst_fileparts(Files{1}), 'data_subspace');
        OutputFile = file_unique(OutputFile);
        bst_save(OutputFile, DataMat, 'v6');
        db_add_data(iStudy, OutputFile, DataMat);
        FilesOut{iOutput} = OutputFile;
    end
end

if Options.isSave
    output = FilesOut;
else
    output = DataMat;
end

end

%% ADDITIONAL FUNCTIONS

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


%% Old function to compute geneig problem

% function varargout = compute_geneig_matlab(CovA, CovB)
% 
% % Solve generalized eigenvalue problem
% [V, D] = eig(CovA, CovB);
% D = diag(D);
% [~, index] = sort(D, 'descend');
% SNR = D(index);
% varargout{2} = V(:, index);
% % invert the filters to obtain the forward patterns
% if nargout>2
%     varargout{3} = pinv(varargout{2})';
% end
% % Other output
% varargout{1} = SNR;
% 
% end

%% Old auto-/cross-covariance function

% function [Cxy, Cxx] = compute_auto_cross_covariance_old(X, y)
% 
% [nChan, nSample, nTrial] = size(X);
% 
% % Normalize data to avoid numerical issues
% X = X ./ mean(mean(std(X,[],2)));
% y = y ./ mean(mean(std(y,[],2)));
% 
% % Initialize covariance matrices
% Cxy = zeros(nChan);
% Cxx = Cxy;
% 
% for iTrial = 1:nTrial
%     
%     % Initialize cross- and auto-correlation
%     xy = zeros(nChan, nSample*2+1);
%     xx = zeros(nChan, nSample*2+1);
%     
%     % compute cross-correlation for each channel, positive lags mean B
%     % following A (e.g. data following reference signal)
%     for iChan = 1:nChan
%         xx(iChan,:) = xcorr(X(iChan,:,iTrial), nSample);
%         xy(iChan,:) = xcorr(X(iChan,:,iTrial), y(:,:,iTrial), nSample);
%     end
%     yy = xcorr(y(:,:,iTrial), nSample);
%     
%     % Ensure zero mean
%     xx = bsxfun(@minus, xx, mean(xx,2));
%     xy = bsxfun(@minus, xy, mean(xy,2));
%     yy = yy - mean(yy);
%     
%     % % Windowing function
%     % win = bst_window('tukey',size(xy,2))';
%     % xy = bsxfun(@times,xy,win);
%     % xx = bsxfun(@times,xx,win);
%     % yy = bsxfun(@times,yy,win);
%     
%     % Compute covariance of auto/cross-correlation
%     Cxy = Cxy + xy * xy';
%     Cxx = Cxx + xx * xx';
%     
% end
% 
% Cxy = Cxy ./ (nSample*2+1)*nTrial;
% Cxx = Cxx ./ (nSample*2+1)*nTrial;
% 
% end

%% Older statistics part in the main function

% % Statistics to determine number of dimensions
% if (Options.isStat)&&(strcmp(pipeline, 'stim_base'))
%     fprintf('Statistics\n');
%     Options.nPerm = 599;
%     Options.nDim = 1;
%     Options.stat0 = 1;
%     Options.func = 'geneig_func';
%     Options.verbose = 1;
%     for iTrial = 1:nTrial
%         inputs1{iTrial} = C1_all(:,:,iTrial);
%         inputs2{iTrial} = C2_all(:,:,iTrial);
%     end
%     pval = permtest_dep2(inputs1, inputs2, Options);
% end
% if (Options.isStat)&&~(strcmp(pipeline, 'stim_base'))
%     if ~(strcmp(pipeline, 'stim_base'))
%         balpoint = sum(svd(C1))/sum(svd(C2));
%     else
%         balpoint = 1;
%     end
%     for iSegment = 1:size(C1_all,3)
%         val1 = diag(Results{2}'*C1_all(:,:,iSegment)*Results{2});
%         val2 = diag(Results{2}'*C2_all(:,:,iSegment)*Results{2});
%         SNR_all(:,iSegment) = val1 ./ val2;
%     end
%     nBoot = 5999;
%     nSeg = size(C1_all,3);
%     clear SNR_boot
%     for iBoot = 1:nBoot
%         bootsample = randsample(1:nSeg,nSeg,1);
%         SNR_boot(:,iBoot) = mean(SNR_all(:,bootsample),2);
%     end
% 
%     pval = 1-sum(SNR_boot>balpoint,2)'/nBoot;
%     iCut = find(diff(pval<.001),1);
%     bootsort = sort(SNR_boot,2);
%     clear ci_spec
%     ci_spec(:,1) = bootsort(:,round(nBoot*.9995));
%     ci_spec(:,2) = bootsort(:,round(nBoot*.0005));    
% end
% if 0
%     if Options.isStat
%         inputs1 = {}; inputs2 = {};
%         for iSegment = 1:size(SNR_all,2)
%             inputs1{iSegment} = SNR_all(:,iSegment);
%             inputs2{iSegment} = 1./SNR_all(:,iSegment);
%         end
%         Options.nPerm = 599;
%         Options.nDim = 1;
%         Options.stat0 = 1;
%         Options.func = inline('mean([input{1}{:}],2)','input','para');
%         para.verbose = 1;
%         pval = permtest_dep2(inputs1, inputs2, Options);
%     end
% end
% if 0
%     if (Options.isStat)&&(strcmp(pipeline, 'freq1_freq2'))
%         Options.flims = [min(Options.PassBand(1),Options.NoiseBand(1)) ...
%             max(Options.PassBand(2),Options.NoiseBand(2))];
%         X = get_data(Files, iGoodChan, Options.Stimulus, [], [], [], []);
%         nSignals = size(Results{2},2);
%         nTrials = size(X,3);
%         for iSignal = 1:nSignals
%             iSignal
%             for iTrial = 1:nTrials
%                 xs = Results{2}(:,iSignal)' * X(:,:,iTrial);
%                 [psdest_tmp, f] = bst_psd(xs, fs, 2, .5);
%                 psdest(iTrial,:) = squeeze(psdest_tmp);
%             end
%             hold on; plot(f(findex),mean(psdest(:,findex)))
%             [pval1, pval2(iSignal), pval3(iSignal)] = boottest_oneoverf(psdest, f, Options.flims, Options.PassBand, 0);
%         end
%     end
% end