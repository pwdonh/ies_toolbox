
SubjectNames = {'RSpert001', 'RSpert002', 'RSpert003', 'RSpert004', ...
                'RSpert005', 'RSpert006', 'RSpert007', 'RSpert008', ...
                'RSpert009', 'RSpert010', 'RSpert011', 'RSpert012', ...
                'RSpert013', 'RSpert014', 'RSpert015', 'RSpert016', ...
                'RSpert017'};
BandPass = {[4 8], [8 13], [13 30]};
noisebp = [2, 100];
nSubj = numel(SubjectNames);
epochtime = [1 3];

%% Compute signal subspace for all subjects

for iSubj = 1:nSubj
    
    EpochCoregFiles = bst_process('CallProcess', 'process_select_files_data', [], [], ...
        'subjectname',   SubjectNames{iSubj}, 'tag', 'start');
    indCond = strcmp('rspert', {EpochCoregFiles.Condition});
    EpochCoregFiles = EpochCoregFiles(indCond);
    
    for iBand = 1:numel(BandPass)
        
        % Process: Compute Signal Subspace: Narrowband
        subspFiles(iSubj,iBand) = bst_process('CallProcess', 'process_subspace_narrowband', EpochCoregFiles, [], ...
            'bandpass',    BandPass{iBand}, ...
            'noisebp',     noisebp, ...
            'stim',        epochtime, ...
            'reg',         0.05, ...
            'split',       0, ...
            'win_length',  [], ...
            'win_overlap', [], ...
            'wfcn',        [], ...
            'isoutf',      0, ...
            'issavecov',   0);            

        % Process: Compute Signal Subspace: Statistics
        bst_process('CallProcess', 'process_subspace_statistics', subspFiles(iSubj,iBand), [], ...
            'niter', 6000);
        
        % Process: Add tag
        bst_process('CallProcess', 'process_add_tag', subspFiles(iSubj,iBand), [], ...
            'tag',    num2str(BandPass{iBand}), ...
            'output', 1);  % Add to comment        
        
        % Process: Scan signal subspace (fast)
        scanFiles(iSubj,iBand) = bst_process('CallProcess', 'process_scan_subspace', subspFiles(iSubj,iBand), [], ...
            'rsspace', 2, ...  % p-value threshold
            'ndim',    [], ...
            'alpha',   0.001, ...
            'snr',     [], ...
            'order',   1, ...  % first
            'whiten',  2, ...  % Noise covariance
            'niter',   1);
        
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
        
        %%% Subspace scanning
        
        % Select subjects who have scan file
        scanFilesOrder = {scanFiles(:,iBand).FileName};
        sigindex = ~cellfun('isempty',scanFilesOrder);
        scanFilesOrder = scanFilesOrder(sigindex);
        
        % Process: Project on default anatomy: volume
        scanFilesGroup = bst_process('CallProcess', 'process_project_sources', scanFilesOrder, [], ...
            'headmodeltype', 'volume');  % MRI volume
        
%         % Process: Delete selected files
%         bst_process('CallProcess', 'process_delete', scanFilesOrder, [], ...
%             'target', 1);  % Delete selected files
        
        % Process: Subspace scanning group statistics
        scanFileAvg = bst_process('CallProcess', 'process_scan_subspace_statistics', scanFilesGroup, [], ...
            'nperm', 6000, ...
            'isinterp', 0, ...
            'niter', 1);
        
%         [~,iStudy] = bst_get('StudyWithCondition', 'Group_analysis/rspert/');
%         avgScanResults{iBand,iOrder} = db_add(iStudy, scanFileAvg.FileName, 1);
        avgScanResults{iBand} = scanFileAvg.FileName;
        
        % Process: Add tag:
        bst_process('CallProcess', 'process_add_tag', avgScanResults{iBand}, [], ...
            'tag',    [num2str(BandPass{iBand}) ' | scanned'], ...
            'output', 1);  % Add to comment
        
        
        %%% Imaging of ratios
        
        % Select subjects who have scan file
        imageFilesSel = {imageFiles(:,iBand).FileName};
        imageFilesSel = imageFilesSel(sigindex);
        
        % Process: Project on default anatomy: volume
        imageFilesGroup = bst_process('CallProcess', 'process_project_sources', imageFilesSel, [], ...
            'headmodeltype', 'volume');  % MRI volume
        
        % Process: Average: Everything
        avgImageResults{iBand} = bst_process('CallProcess', 'process_average', imageFilesGroup, [], ...
            'avgtype',         1, ...  % Everything
            'avg_func',        1, ...  % Arithmetic average:  mean(x)
            'weighted',        0, ...
            'scalenormalized', 0);        
        
        % Process: Add tag:
        bst_process('CallProcess', 'process_add_tag', avgImageResults{iBand}, [], ...
            'tag',    [num2str(BandPass{iBand}) ' | imaged'], ...
            'output', 1);  % Add to comment
        
end


