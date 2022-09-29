function Create_ECoG_Epoch_RD()

    % RepSup Subjects 
    sub_str = {'357', '362', '369', '372', '376'};
    sub_str = {'357'};
    session2_str = {'', '', '', '', '', ''};
    
    for i=1:length(sub_str)
% % %         DoWork(sub_str{i}, session2_str{i}, 'stimuli_1');
% % %         DoWork(sub_str{i}, session2_str{i}, 'onset_1');
        DoWork(sub_str{i}, session2_str{i}, 'stimuli_12');
        DoWork(sub_str{i}, session2_str{i}, 'onset_12');

    end

end

function outlier_bounds = DoWork(sub_str, session2_str, data_type)
    plot_data = false;
    remove_outliers_flag = true;
    save_flag = false;
    exclude_unused = false;
    only_active = false; % Active designation is now done in SubjectClass
    baseline_from_1 = true;


    % Variables:

    % RepSup Subjects 
%     sub_str = '357'; % {'352', '357', '362', '369', '372', '376'}
%     session2_str = '';
    local_onset_flag = true;

    freqBand = 'highgamma'; % {'highgamma', 'broadband', 'beta'}
    if strcmp(freqBand, 'highgamma')
        band_short = 'HG';
    elseif strcmp(freqBand, 'beta')
        band_short = 'BT';
    else
        error('Implement');
    end

    preprocessing_method = 'Hilbert'; % {'Hilbert', 'MultiTaper', 'Filter', 'Wavelet'}

%     condition = 'stimuli_1'; % {'stimuli_1', 'stimuli_2', 'stimuli_12', ...
%                                 %  'onset_1', 'onset_2', 'onset_12'}


    switch data_type
        case 'stimuli_1'
            pre_time = 1;
            post_time = 2;
            baseline_time = 0.5;
        case 'onset_1'
            pre_time = 2;
            post_time = 1;
            baseline_time = 0.5;
        %?? Added these last two cases
        case 'stimuli_12'
            pre_time = 1;
            post_time = 2;
            baseline_time = 0.5;
        case 'onset_12'
            pre_time = 2;
            post_time = 1;
            baseline_time = 0.5;
    end

    if local_onset_flag
        local_str = 'Local';
    else
        local_str = '';
    end

    if post_time == 2
        add_str = ''; % additional string to add to end of file names
    else
        add_str = 'postNot2s';
    end
    if baseline_from_1
        add_str = [add_str '_baseline1'];
    else
        add_str = [add_str '_baseline2'];
    end


% % %     base_path = error('Point to data directory here');
% %     base_path = 'R:\projectnb\busplab\Experiments\ECoG_Preprocessed';
    base_path = '/projectnb/busplab/Experiments/ECoG_Preprocessed_RD';
    if ~save_flag
        fprintf('Not saving results.\n');
    else
        if baseline_from_1
            save_path = fullfile(base_path, 'LocalEpoched', 'BSL1', ['S' sub_str]);
        else
            save_path = fullfile(base_path, 'LocalEpoched', 'BSL12', ['S' sub_str]);
        end
        if ~exist(save_path,'dir')
            mkdir(save_path);
            fprintf('Creating: %s...\n',save_path);
        end
        fprintf('Saving to: %s...\n', save_path);
        fig_save_path = fullfile(save_path, ['Figures_' datestr(now, 'yyyymmdd')]);

        % Save location
        if ~exist(fig_save_path, 'dir') && plot_data
            fprintf('Creating %s...\n', fig_save_path);
            mkdir(fig_save_path);
        end
    end

    % Add paths

    % addpath(genpath('EEG-Clean-Tools-master'));
    % addpath(genpath('eeglab13_5_4b'));
% % %     addpath('PreprocessTools\');
% % %     addpath(fullfile(base_path, 'RawData'));


    % Load files

    % Load data file
    load(fullfile(base_path, 'LocalProcessed', ['S' sub_str], ['ECoGALL' session2_str '_' preprocessing_method '_' band_short '.mat']));

    % Load event file
    load(fullfile(base_path, 'LocalProcessed', ['S' sub_str], [local_str 'OnsetTable' session2_str '.mat']));

    % Load chan_list
    load(fullfile('/project/busplab/software/ecog/data/LocalEpoched_chan_ids.mat'));
    chan_list = chan_list.(['S' sub_str]);
        
    if exclude_unused
        chans = Unused_Electrodes(str2double(sub_str));       
        chan_list = chan_list(~ismember(chan_list, chans));
    end
    
    
% % %     % Subject specific files
% % %     sub_vars = subject_vars(sub_str);
% % %     % chan_list = sub_vars.valid_electrodes;
% % %     chan_list = setdiff(sub_vars.valid_electrodes, sub_vars.bad_chans);


    % Variables

    num_sessions = length(ECoGALL);

    sessions_to_process = 1:num_sessions;

    epoch_data_raw = [];
    bl_epoch_data_raw = [];
    event_info_raw = [];
    bl_event_info_raw = [];
    ids_no_nan = {};
    bad_channels = cell(1, num_sessions);
    for sess = sessions_to_process

        switch preprocessing_method
            case 'Hilbert'
                X = ECoGALL{sess}.data_hilbert;
            case 'Filter'
                X = ECoGALL{sess}.data_filt;
            case 'MultiTaper'
                X = ECoGALL{sess}.data_MT;
            otherwise
                X = ECoGALL{sess}.data;
        end
        onsetTable = OnsetTable{sess};

        % Channels to evaluate:
        fprintf('Additional Bad channel removal here?\n');        

    
    
        X = X(chan_list, :);

    %     EEG.Fs = EEG.Fs(good_chans);
    %     EEG.data = EEG.data(good_chans, :);
    %     EEG.nchan = length(good_chans);
    %     EEG.data_hilbert = EEG.data_hilbert(good_chans, :);
    %     EEG.badChan = intersect(good_chans, EEG.badChan);

        num_chans = size(X, 1);
        num_samples = size(X, 2);

        if sess > 1
            if Fs ~= ECoGALL{sess}.srate
                error('Sampling rates should be the same');
            end
        else
            Fs = ECoGALL{sess}.srate;

            dt = 1/Fs;

            pre_samples = pre_time * Fs;
            post_samples = post_time * Fs; % Includes onset sample

            % Defining baseline period
            baseline_samples = baseline_time * Fs;
            baseline_period = 1:baseline_samples;

        end



        % Create Epochs

        % Extract Events (gets sample number for event)
        % Columns are setup as follows (6 columns):
        % Presentation of 1st word, presentation of 2nd word, 
        % speech onset of 1st word, speech onset of 2nd word, 
        % speech offset of 1st word, speech offset of 2nd word

        [events, sess_event_info, ids, bl_events, bl_sess_event_info] = GetEpochTimesRD(onsetTable, data_type, baseline_from_1);
        
        ids_no_nan{end+1} = ids;
        
        num_events = length(events);
        extractInds = -pre_samples : post_samples-1;
        extractTime = extractInds/Fs;

        epoch = nan(num_chans, pre_samples+post_samples, num_events);
%         if baseline_from_1
            bl_epoch = nan(num_chans, length(baseline_period), num_events);
            assert(mod(num_events,length(bl_events)) == 0);
%         else
            
%         end
        for event = 1:num_events

            inds2use = events(event) + extractInds;
            if strcmp(data_type, 'stimuli_12') 
%                 bl_inds2use = bl_events(event)  - pre_samples + baseline_period - 1;
                bl_inds2use = bl_events(event) - fliplr(baseline_period);
            elseif strcmp(data_type, 'onset_12')
                if baseline_from_1
                    bl_inds2use = onsetTable(ids(round(event/2)), 1) - fliplr(baseline_period);
                else
                    try
                        bl_inds2use = onsetTable(ids(round(event/2)), mod(event-1,2) + 1) - fliplr(baseline_period);
                    catch
                        pause
                    end
                end
            else
                fprintf('New Condition\n');
                pause;
            end

            epoch(:, :, event) = X(:, inds2use);
%             if baseline_from_1
            try
                bl_epoch(:, :, event) = X(:, bl_inds2use);
            catch
                pause
            end
%             end

        end

        epoch_data_raw = cat(3, epoch_data_raw, epoch);
        event_info_raw = [event_info_raw; sess_event_info];

%         if 
            bl_epoch_data_raw = cat(3, bl_epoch_data_raw, bl_epoch);
            bl_event_info_raw = [bl_event_info_raw; bl_sess_event_info];
%         end


    end




    % Create the preprocessor class
    preprocessor = PreprocessClass(epoch_data_raw, chan_list, Fs);







    % calculate zscore

%     epoch_data_z = preprocessor.z_score(epoch_data_raw, baseline_period);
%     if baseline_from_1
    epoch_data_z = preprocessor.z_score_with_bl_data(epoch_data_raw, bl_epoch_data_raw);
%     end


    % check for outliers (protocol from Megan's filter_ECoG_trials_across)
    if remove_outliers_flag
        num_outliers = zeros(size(epoch_data_z,1),1);
        for i_ch = 1:size(epoch_data_z,1)
            % Find peaks and median
            ch = epoch_data_z(i_ch,:,:);
            ch_peaks = squeeze(max(ch,[],2))';
            ch_troughs = squeeze(min(ch,[],2))';
            
            med_ch = median(ch,3);
            med_peak = median(ch_peaks);        
            med_trough = median(ch_troughs);

            % Scale difference in peaks from median
    %         ch_peaks_scaled = (ch_peaks - med_peak) ./ med_peak;
    %         ch_scaled = (ch - med_ch) ./ med_ch;

            % Calculate threshold
            lb = 4*std(rmoutliers(ch_troughs, 'percentiles', [10 90]));
            ub = 4*std(rmoutliers(ch_peaks, 'percentiles', [10 90]));


            % Identify Outliers
            outlier_bounds{i_ch} = [med_trough - lb, med_peak + ub];
            idx_outliers = find(ch_peaks >= med_peak + ub | ch_troughs <= med_trough - lb);
            idx_not_outliers = 1:size(epoch_data_z,3);
            idx_not_outliers(idx_outliers) = [];
            fprintf('%d < x < %d\n', med_trough - lb, med_peak + ub)
            num_outliers(i_ch) = size(idx_outliers,2);
            for i = 1:size(idx_outliers,2)
               fprintf('Outlier %3.d/%3.d: Ch - %d\n',i, size(idx_outliers,2), i_ch)

               % Accept/Reject Outliers
%                close; figure; hold on;
               if ~any(isnan(epoch_data_z(i_ch,:,idx_outliers(i))))
%                    plot(mean(epoch_data_z(i_ch,:,idx_not_outliers),3,'omitnan'))
%                    plot(epoch_data_z(i_ch,:,idx_outliers(i)))
%                    legend('Mean Trace for Channel', 'Outlier Trace')

    %                str_in = input('Keep trial in? (y/n): ', 's');
    %                if ~strcmpi(str_in,'y')
                       if mod(idx_outliers(i),2)
                           epoch_data_z(i_ch,:,idx_outliers(i)) = nan(size(epoch_data_z(i_ch,:,idx_outliers(i))));
                           epoch_data_z(i_ch,:,idx_outliers(i)+1) = nan(size(epoch_data_z(i_ch,:,idx_outliers(i))));
                       else
                           epoch_data_z(i_ch,:,idx_outliers(i)) = nan(size(epoch_data_z(i_ch,:,idx_outliers(i))));
                           epoch_data_z(i_ch,:,idx_outliers(i)-1) = nan(size(epoch_data_z(i_ch,:,idx_outliers(i))));
                       end

    %                end
               end

            end



        end
    end
    
    % 
    if only_active
        [chan_list, data_idx] = Active_Electrodes(str2num(sub_str), data_type, epoch_data_z, chan_list);
        epoch_data_z = epoch_data_z(data_idx, :, :);
    end

    preprocessed_data = AnalysisClass(epoch_data_z, chan_list, Fs);
    preprocessed_data.event_sample = pre_samples;
    electrodes2lookAt = 1:size(epoch_data_raw, 1);

    %

    if plot_data

        addpath(genpath('Tools'));

        fig_0 = preprocessed_data.SetupGridPlot(electrodes2lookAt');

        chans = fig_0.plot_ref;
        % TODO: Get plotting in correct order

        max_val = 10; %max(max(max(epoch_data_z)));
        min_val = -10;%min(min(min(epoch_data_z)));
        for ind = 1:numel(fig_0.plot_ref)

            if isnan(chans(ind))
                continue;
            end

            % Get the axis going from top left across first
            aH = fig_0.plot_axes_handles(ind);
            % Get the channel location for same order
            chan_ind = fig_0.plot_ref(ind);

            axes(aH);
            box(aH, 'off');
            hold(aH, 'all');

            imagesc(squeeze(epoch_data_z(chans(ind), :, :))');
            set(gca, 'CLim', [min_val max_val]);
            colormap(flipud(brewermap([], 'RdBu')))
        end


    %

        f1_fig = preprocessed_data.PlotERP();
        f1 = f1_fig.handle;


    end



    %


    if save_flag 

        % Data
        event_info = struct('events', event_info_raw, 'ids_no_nan', {ids_no_nan}, 'Fs', Fs, 'condition', data_type, ...
                            'bl_events', bl_event_info_raw);
        save_file = fullfile(save_path, ['Epoch' session2_str '_' data_type '_' preprocessing_method '_' band_short '.mat']);
        save(save_file, 'preprocessed_data', 'event_info');

        % Save figures

        % ERSP
        if plot_data
            save_filename = ['ERSP_' data_type];
            set(f1, 'Position', [180, 100, 1175, 650])
            fname1 = fullfile(fig_save_path, [save_filename session2_str '.fig']);
            fname2 = fullfile(fig_save_path, [save_filename session2_str '.png']);
            fprintf('Saving %s and %s...\n', fname1, fname2);
    %         saveas(f1, fname1, 'fig');
            saveas(f1, fname2, 'png');
            close(f1);

            % Trial Image
            f0 = fig_0.handle;
            save_filename = ['TrialImage_' data_type];
            set(f0, 'Position', [180, 100, 1175, 650])
            fname1 = fullfile(fig_save_path, [save_filename session2_str '.fig']);
            fname2 = fullfile(fig_save_path, [save_filename session2_str '.png']);
            fprintf('Saving %s and %s...\n', fname1, fname2);
    %         saveas(f0, fname1, 'fig');
            saveas(f0, fname2, 'png');
            close(f0);
        end

    end

end

function chans = Unused_Electrodes(sub_str)
% The listed electrodes are excluded
switch sub_str
    case 357
        chans = [5 6 61:64 71:76 95 104 105 161:224];
    case 362
        chans = [5,6,63,64,93:96,207:224];        
    case 372
        chans = [1:66 85:96 109:128];
    otherwise
        chans = [];    
end
end

function [ae, data_idx] = Active_Electrodes(sub_str, data_type, data, chan_list)
win = 250;
start = 2000;
stop = 3000;
threshold = 1;
nDur = 100;

% Truncate to window
if strcmp(data_type, 'stimuli_12')
    data = data(:,start:stop,:);
else
    data = data(:,start:stop,:); % Change for onset? Would need to find SO time in onset_12
end

% Load conditions
load('/projectnb/busplab/Experiments/ECoG_Preprocessed_RD/filelist.mat', 'filelist');
trials = filelist.trials(filelist.subject == sub_str);
trials = trials{1,1};
conds = sort(unique(trials.pair_comparison));

% Split dataset by conditions and word, then average
for i_cond = 1:size(conds,1)% different, flipped, identical
    idx = strcmp(trials.pair_comparison, conds{i_cond});
    temp = data(:, :, idx);
    
    w1(:, :, i_cond) = mean(temp(:, :, 1:2:end), 3, 'omitnan');
    w2(:, :, i_cond) = mean(temp(:, :, 2:2:end), 3, 'omitnan'); 
end

% Process
chans = cell(size(conds,1),1); data_chans = cell(size(conds,1),1);
temp = [];
for i_cond = 1:size(conds,1)
    % word 1
    for i_win = 1:size(w1,2)-win
        temp(:,i_win) = abs(mean(w1(:,i_win+win,i_cond),2,'omitnan')) > threshold;
    end
    data_chans1 = find(sum(temp,2) > nDur);
    chans1 = chan_list(data_chans1);
    
    
    % word 2
    for i_win = 1:size(w2,2)-win
        temp(:,i_win) = abs(mean(w2(:,i_win+win,i_cond),2,'omitnan')) > threshold;
    end
    data_chans2 = find(sum(temp,2) > nDur);
    chans2 = chan_list(data_chans2);
    
    % channels in either w1 and w2
    data_chans{i_cond} = union(data_chans1, data_chans2);
    chans{i_cond} = union(chans1, chans2);
end

% channels in both w1 and w2 across any condition

data_idx = union(union(data_chans{1}, data_chans{2}), data_chans{3});
ae = union(union(chans{1}, chans{2}), chans{3});
end