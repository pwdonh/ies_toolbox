
%% Analysis parameters

SubjectNames = {'RSpert001', 'RSpert002', 'RSpert003', 'RSpert004', ...
                'RSpert005', 'RSpert006', 'RSpert007', 'RSpert008', ...
                'RSpert009', 'RSpert010', 'RSpert011', 'RSpert012', ...
                'RSpert013', 'RSpert014', 'RSpert015', 'RSpert016', ...
                'RSpert017'};
BandPass = {[50 85], [13 30]};
nSubj = numel(SubjectNames);
stimtime = [1 3];
basetime = [-2 0];

%% Subject-level analysis

for iBand = 1:numel(BandPass)
    
    for iSubj = 1:numel(SubjectNames)
        
        EpochCoregFiles = bst_process('CallProcess', 'process_select_files_data', [], [], ...
            'subjectname',   SubjectNames{iSubj}, 'tag', 'start');
        indCond = strcmp('rspert', {EpochCoregFiles.Condition});
        EpochCoregFiles = EpochCoregFiles(indCond);
        
        % Process: Compute Signal Subspace: Stimulus over Baseline
        subspFiles(iSubj,iBand) = bst_process('CallProcess', 'process_subspace_induced', EpochCoregFiles, [], ...
            'bandpass',    BandPass{iBand}, ...
            'stim',        stimtime, ...
            'base',        basetime, ...
            'ispca',       0, ...
            'noisebp',     [0 0], ...
            'reg',         0.05, ...
            'split',       0, ...
            'win_length',  4, ...
            'win_overlap', 50, ...
            'wfcn',        'tukey', ...
            'isoutf',      0, ...
            'issavecov',   1);

        % Process: Add tag:
        bst_process('CallProcess', 'process_add_tag', subspFiles(iSubj,iBand), [], ...
            'tag',    num2str(BandPass{iBand}), ...
            'output', 1);  % Add to comment        
        
        bst_process('CallProcess', 'process_subspace_statistics', subspFiles(iSubj,iBand), [], ...
            'niter', 600);
        
        % Remove single-trial cov. matrices to save storage space
        sMat = in_bst_data(subspFiles(iSubj,iBand).FileName);
        sMat = rmfield(sMat, 'C1_all');
        sMat = rmfield(sMat, 'C2_all');
        bst_save(file_fullpath(subspFiles(iSubj,iBand).FileName), sMat);
%         
        % Scan signal subspace
        for iOrder = 1:2 % stimulus > baseline or baseline > stimulus
            if iOrder==1
                isSig = sMat.pval(1)<.05;
                nDim = find(diff(cumsum(sMat.pval' < 0.05))==0,1);
            else
                isSig = sMat.pval(end)<.05;
                nDim = find(diff(cumsum(flip(sMat.pval)'<0.05))==0,1);
            end
            if isSig
                % Process: Scan signal subspace (fast)
                scanFiles(iSubj,iOrder,iBand) = bst_process('CallProcess', 'process_scan_subspace', subspFiles(iSubj,iBand), [], ...
                    'rsspace', 2, ...  % p-value threshold
                    'ndim',    [], ...
                    'alpha',   0.05, ...
                    'snr',     [], ...
                    'order',   iOrder, ...  % first
                    'whiten',  2, ...  % Noise covariance
                    'niter',   1);
                
                % This will extract the component time-series:
%                 % Process: Apply spatial filter
%                 bst_process('CallProcess', 'process_apply_spfilter2', subspFiles(iSubj,iBand), EpochCoregFiles, ...
%                     'ncomp',       nDim, ...
%                     'order',       iOrder, ...  % first
%                     'matrix',      2, ...  % Components (one matrix per epoch)
%                     'comment',     [subspFiles(iSubj,iBand).Comment ' | signals'], ...
%                     'isnoisenorm', 0);

            end
        end
        
        % Image signal subspace
        sStudy = bst_get('Study',subspFiles(iSubj,iBand).iStudy);
        isKernel = ~cellfun('isempty', strfind({sStudy.Result.FileName},'KERNEL'));
        isLink = [sStudy.Result.isLink];
        isVolume = strcmp({sStudy.Result.HeadModelType},'volume');
        KernelFile = sStudy.Result(isKernel&(~isLink)&isVolume).FileName;
        resultFile = ['link|' KernelFile '|' subspFiles(iSubj,iBand).FileName];
        
        imageFiles(iSubj,iBand) = bst_process('CallProcess', 'process_image_subspace', resultFile, []);
        
    end
    
end

%% Group analysis

for iBand = 1:numel(BandPass)
    for iOrder = 1:2
        
        %%% Subspace scanning
        
        % Select subjects who have scan file
        scanFilesOrder = {scanFiles(:,iOrder,iBand).FileName};
        sigindex(:,iOrder) = ~cellfun('isempty',scanFilesOrder);
        scanFilesOrder = scanFilesOrder(sigindex(:,iOrder));
        
        % Process: Project on default anatomy: volume
        scanFilesGroup = bst_process('CallProcess', 'process_project_sources', scanFilesOrder, [], ...
            'headmodeltype', 'volume');  % MRI volume
        
%         % Process: Delete selected files
%         bst_process('CallProcess', 'process_delete', scanFilesOrder, [], ...
%             'target', 1);  % Delete selected files
        
        % Process: Subspace scanning group statistics
        scanFileAvg = bst_process('CallProcess', 'process_scan_subspace_statistics', scanFilesGroup, []);
        
%         [~,iStudy] = bst_get('StudyWithCondition', 'Group_analysis/rspert/');
%         avgScanResults{iBand,iOrder} = db_add(iStudy, scanFileAvg.FileName, 1);
        avgScanResults{iBand,iOrder} = scanFileAvg.FileName;
        
        % Process: Add tag:
        bst_process('CallProcess', 'process_add_tag', avgScanResults{iBand,iOrder}, [], ...
            'tag',    [num2str(BandPass{iBand}) ' | order ' num2str(iOrder) ' | scanned'], ...
            'output', 1);  % Add to comment
        
        
        %%% Imaging of ratios
        
        % Select subjects who have scan file
        imageFilesSel = {imageFiles(:,iBand).FileName};
        imageFilesSel = imageFilesSel(sigindex(:,iOrder));
        
        % Process: Project on default anatomy: volume
        imageFilesGroup = bst_process('CallProcess', 'process_project_sources', imageFilesSel, [], ...
            'headmodeltype', 'volume');  % MRI volume
        
        % Process: Subspace scanning group statistics
        imageFileAvg = bst_process('CallProcess', 'process_image_subspace_statistics', imageFilesGroup, []);
        
%         [~,iStudy] = bst_get('StudyWithCondition', 'Group_analysis/method_paper_induced/');
%         avgImageResults{iBand,iOrder} = db_add(iStudy, imageFileAvg.FileName, 1);
        avgImageResults{iBand,iOrder} = imageFileAvg.FileName;
        
        % Process: Add tag:
        bst_process('CallProcess', 'process_add_tag', avgImageResults{iBand,iOrder}, [], ...
            'tag',    [num2str(BandPass{iBand}) ' | order ' num2str(iOrder) ' | imaged'], ...
            'output', 1);  % Add to comment
        
    end
end







