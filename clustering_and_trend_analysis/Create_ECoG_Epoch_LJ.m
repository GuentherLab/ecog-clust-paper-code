function Create_ECoG_Epoch()

    % RepSup Subjects 
    sub_str = {'352', '357', '362', '369', '372', '376'};
    session2_str = {'', '', '', '', '', ''};
    
    for i=3%:length(sub_str)
        DoWork(sub_str{i}, session2_str{i}, 'stimuli_1');
        DoWork(sub_str{i}, session2_str{i}, 'onset_1');
    end

end

function DoWork(sub_str, session2_str, condition)

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

    plot_data = true;
    save_flag = true;

    baseline_from_1 = false;
    switch condition
        case 'stimuli_1'
            pre_time = 1;
            post_time = 2;
            baseline_time = 0.5;
        case 'onset_1'
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
    end


%     base_path = error('Point to data directory here');
    base_path = '/projectnb/busplab/Experiments/ECoG_Preprocessed';
    
    
    if ~save_flag
        fprintf('Not saving results.\n');
    else
        save_path = fullfile(base_path, 'LocalEpoched', ['S' sub_str]);
        fprintf('Saving to: %s...\n', save_path);
        fig_save_path = fullfile(base_path, 'LocalEpoched', ['S' sub_str], ['Figures_' datestr(now, 'yyyymmdd')]);

        % Save location
        if ~exist(fig_save_path, 'dir') && plot_data
            fprintf('Creating %s...\n', fig_save_path);
            mkdir(fig_save_path);
        end
    end

    % Add paths

    % addpath(genpath('EEG-Clean-Tools-master'));
    % addpath(genpath('eeglab13_5_4b'));
    addpath('PreprocessTools\');
    addpath(fullfile(base_path, 'RawData'));


    % Load files

    % Load data file
    load(fullfile(base_path, 'LocalProcessed', ['S' sub_str], ['ECoGALL' session2_str '_' preprocessing_method '_' band_short '.mat']));

    % Load event file
    load(fullfile(base_path, 'LocalProcessed', ['S' sub_str], [local_str 'OnsetTable' session2_str '.mat']));


    % Subject specific files
    sub_vars = subject_vars(sub_str);
    % chan_list = sub_vars.valid_electrodes;
    chan_list = setdiff(sub_vars.valid_electrodes, sub_vars.bad_chans);


    % Variables

    num_sessions = length(ECoGALL);

    sessions_to_process = 1:num_sessions;

    epoch_data_raw = [];
    bl_epoch_data_raw = [];
    event_info_raw = [];
    bl_event_info_raw = [];
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

        [events, sess_event_info, bl_events, bl_sess_event_info] = GetEpochTimes(onsetTable, condition, baseline_from_1);

        num_events = length(events);
        extractInds = -pre_samples : post_samples-1;
        extractTime = extractInds/Fs;

        epoch = nan(num_chans, pre_samples+post_samples, num_events);
        if baseline_from_1
            bl_epoch = nan(num_chans, length(baseline_period), num_events);
            assert(num_events == length(bl_events));
        end
        for event = 1:num_events

            inds2use = events(event) + extractInds;
            if baseline_from_1
                bl_inds2use = bl_events(event)  - pre_samples + baseline_period - 1;
            end

            epoch(:, :, event) = X(:, inds2use);
            if baseline_from_1
                bl_epoch(:, :, event) = X(:, bl_inds2use);
            end

        end

        epoch_data_raw = cat(3, epoch_data_raw, epoch);
        event_info_raw = [event_info_raw; sess_event_info];

        if baseline_from_1
            bl_epoch_data_raw = cat(3, bl_epoch_data_raw, bl_epoch);
            bl_event_info_raw = [bl_event_info_raw; bl_sess_event_info];
        end


    end




    % Create the preprocessor class
    preprocessor = PreprocessClass(epoch_data_raw, chan_list, Fs);







    % calculate zscore

    epoch_data_z = preprocessor.z_score(epoch_data_raw, baseline_period);
    if baseline_from_1
        epoch_data_z = preprocessor.z_score_with_bl_data(epoch_data_raw, bl_epoch_data_raw);
    end


    % 

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
        event_info = struct('events', event_info_raw, 'Fs', Fs, 'condition', condition, ...
                            'bl_events', bl_event_info_raw);
        save_file = fullfile(save_path, ['Epoch' session2_str '_' condition '_' preprocessing_method '_' band_short '.mat']);
        save(save_file, 'preprocessed_data', 'event_info');

        % Save figures

        % ERSP
        save_filename = ['ERSP_' condition];
        set(f1, 'Position', [180, 100, 1175, 650])
        fname1 = fullfile(fig_save_path, [save_filename session2_str '.fig']);
        fname2 = fullfile(fig_save_path, [save_filename session2_str '.png']);
        fprintf('Saving %s and %s...\n', fname1, fname2);
%         saveas(f1, fname1, 'fig');
        saveas(f1, fname2, 'png');
        close(f1);

        % Trial Image
        f0 = fig_0.handle;
        save_filename = ['TrialImage_' condition];
        set(f0, 'Position', [180, 100, 1175, 650])
        fname1 = fullfile(fig_save_path, [save_filename session2_str '.fig']);
        fname2 = fullfile(fig_save_path, [save_filename session2_str '.png']);
        fprintf('Saving %s and %s...\n', fname1, fname2);
%         saveas(f0, fname1, 'fig');
        saveas(f0, fname2, 'png');
        close(f0);

    end

end