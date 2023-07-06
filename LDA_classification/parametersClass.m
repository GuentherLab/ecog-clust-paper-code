classdef parametersClass < dynamicprops
    % parametersClass defines container of relevant info for analysis parameters
    %   Pass: (ecog_preprocessed_path, event_duration, window, stride,
    %   grouping_var, topN_feat_pool, topN_feat_indiv).
    %   Outputs an object with params accessible through dot indexing.

    properties (SetAccess = protected)
        ecog_preprocessed_path % Path to Preprocessed ECoG data
        epoched_data_path % Immutable, defined based on expected file structure of parent directory
        ml_path % Immutable, defined based on expected file structure of parent directory
        algo_data_path % Immutable, defined based on expected file structure of parent directory
        LDA_path % Immutable, defined based on expected file structure of parent directory
        figs_quickref_path % Immutable, created during pipeline for easy access to figures generated
        times_path % Generated path based on parameters passed to constructor
        times_label % Shortened version of times_path used for IDing some generated files
        grouping_path % Generated path to contain results spectific to a grouping parameter
        grouping_label % Shortened version of grouping_path used for IDing some generated files

        full_run_label % Combined label of the form: (times_label)_(grouping_label)
        electrodes_table % Electrodes table imported from source electrodes.mat file
        grouping_variable % User-defined (but immutable after intialization) variable of electrodes_table to use as a grouping parameter. 'none' also valid

        timing_mode % dynamic or strict (or Not specified, default to strict and assert numeric values were passed in timing_parameters struct)
        fixed_number_features % if timing mode is dynamic, this is included to dictate the number of features required/expected based on timings
        topN_feat_pool % Top N features to select from the mRMR results in the pooled data analyses
        topN_feat_indiv % Top N features to select from the mRMR results in the individual data analyses
    end

    properties
        event_duration % The event duration specifies how much of the epoched signal will be used in the discretization process
        window % The window size (ms) of the moving average filter
        stride % The stride length (ms) of the MVA window
        dynamic_window_timings % Table containing cluster-specific timings (event, window, stride)
        current_group_value % Current grouping variable value, changes per iteration (eg. if grouping_variable is depth_type, current_group_value will cycle through 'depth', 'surface')
        pooled_electrode_feature_set % The pooled electrode feature set associated with the current_group_value
        individual_electrode_feature_set % The individual electrode feature sets associated with the current_group_value
        auc_results_path % Generated path for saving Permutation Test results
    end

    properties (Constant)
        p_value = 0.05; % p-value for determining significance
        number_of_nodes = 50; % number of SCC nodes to use in permutation test
        class_names = {'word_name', 'consonants_name', 'vowel_name'}; % These correspond to the groups of labels
        class_names_formal = {'Word Name', 'Consonants Name', 'Vowel Name'}; % "Formal" labels for titling figures, legends, etc
        sub_nums = [357, 362, 369, 372, 376]; % The subject identifiers as integers
        sub_nums_formal = {'S357', 'S362', 'S369', 'S372', 'S376'}; % The subject identifiers as strings
        sampling_frequency = 1000; % The data was captured at 1kHz
    end
    methods
        % Constructor:
        function params = parametersClass(ecog_preprocessed_path, ...
                timing_parameters, grouping_var, ...
                topN_feat_pool, topN_feat_indiv)
            % Pass: (ecog_preprocessed_path, event_duration, window, stride,
            %     grouping_var, topN_feat_pool, topN_feat_indiv).
            if nargin > 0
                params.check_contents(timing_parameters);
                params.ecog_preprocessed_path = ecog_preprocessed_path;
                params.epoched_data_path = fullfile(ecog_preprocessed_path, 'LocalEpoched');
                params.ml_path = fullfile(ecog_preprocessed_path, 'MachineLearning');
                params.algo_data_path = fullfile(params.ml_path, 'AlgoData');
                params.LDA_path = fullfile(params.ml_path, 'LDA');
                params.figs_quickref_path = fullfile(params.LDA_path, 'figs_quickref');
                params.times_path = fullfile(params.LDA_path, params.times_label);
                params.grouping_label = sprintf('grouped_by_%s', grouping_var);
                params.grouping_path = fullfile(params.times_path, params.grouping_label);
                params.full_run_label = strcat(params.times_label, '_', params.grouping_label);
                params.grouping_variable = grouping_var;
                params.topN_feat_pool = topN_feat_pool;
                params.topN_feat_indiv = topN_feat_indiv;
                params.dynamic_window_timings = readtable('data/dynamic_window_per_cluster.xlsx');
                params.electrodes_table = load('electrodes.mat').electrodes;
            end

        end

        function obj = check_contents(obj, s)

            if isfield(s, 'timing_mode') && strcmp(s.timing_mode, 'dynamic')
                obj.timing_mode = 'dynamic';
                assert(isfield(s, 'fixed_number_features'), 'You must supply a desired number of features to create across cluster timings.');
                disp('Dynamic mode selected! Ignoring values for event_duration, window, and stride and using timings in dynamic_window_per_cluster.xlsx.')
                obj.event_duration = NaN;
                obj.window = NaN;
                obj.stride = NaN;
                obj.fixed_number_features = s.fixed_number_features;
                obj.times_label = 'dynamic';
            else
                obj.timing_mode = 'strict';
                assert(all(ismember({'event_duration', 'window', 'stride'}, fieldnames(s))) && all([~isnan(s.event_duration), ~isnan(s.window), ~isnan(s.stride)]), 'Strict mode selected. You must supply event_duration, window, and stride.');
                disp('Using strict windowing!')
                obj.event_duration = s.event_duration;
                obj.window = s.window;
                obj.stride = s.stride;
                obj.fixed_number_features = s.event_duration / s.stride - 1;
                obj.times_label = sprintf('e%d_w%d_s%d', s.event_duration, s.window, s.stride);
            end

        end
        function [database, electrodes_database] = generate_database(obj)
            % generate_database generates database of subject level data, database of electrode level data.
            %    p = parametersClass(args)
            %    [db, edb] = p.generate_database
            %    db is a table:
            %
            %                 sub_num    data_per_electrode    onset_data_per_electrode
            %                 _______    _______________________    ________________________
            %
            %                   357           {153×1 cell}                {153×1 cell}
            %                   362           {198×1 cell}                {198×1 cell}
            %                   369           {250×1 cell}                {250×1 cell}
            %                   372           {214×1 cell}                {214×1 cell}
            %                   376           {221×1 cell}                {221×1 cell}
            %
            %    edb has the same structure as electrodes.mat (with a few additional columns)
            %    More importantly, edb has a feat_set for each electrode,
            %    based on values in parametersClass.

            if ~isfile(fullfile(obj.times_path, 'database.mat'))
                database = table();

                for sub_num_formal_idx = 1:length(obj.sub_nums_formal)
                    sub_num = obj.sub_nums(sub_num_formal_idx);
                    sub_num_formal = obj.sub_nums_formal{sub_num_formal_idx};

                    raw_data_filelist = dir([fullfile(obj.epoched_data_path, sub_num_formal) '*/*.mat']);
                    raw_data_files = {raw_data_filelist.name};

                    % Don't remember the importance of this block - LJ
                    % BEGIN DEPRECATED  - VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
                    % stim_electrodes_table = obj.electrodes_table(obj.electrodes_table.subject == sub_num & obj.electrodes_table.type == 2, :);
                    % stim_electrodes_list = unique(stim_electrodes_table.idx_LocalProcessed);
                    %                     stim_electrodes_table = stim_electrodes_table(stim_electrodes_list, :);
                    % onset_electrodes_table = obj.electrodes_table(obj.electrodes_table.subject == sub_num & obj.electrodes_table.type == 1, :);
                    % onset_electrodes_list = unique(onset_electrodes_table.idx_LocalProcessed);
                    %                     onset_electrodes_table = onset_electrodes_table(onset_electrodes_list, :);
                    % END DEPRECATED    - ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

                    database_row_temp = table();
                    database_row_temp.('sub_num') = sub_num;
                    for file_idx = 1:numel(raw_data_files)

                        if contains(raw_data_filelist(file_idx).name, 'stim')
                            alignment_type = 'stim';
                        elseif contains(raw_data_filelist(file_idx).name, 'onset')
                            alignment_type = 'onset';
                        else
                            continue;
                        end
                        table_col = sprintf('%s_data_per_electrode', alignment_type);
                        analysis_class = load(fullfile(raw_data_filelist(file_idx).folder, raw_data_filelist(file_idx).name), 'preprocessed_data');

                        preprocessed_data = analysis_class.preprocessed_data.data;
                        if strcmp(obj.timing_mode, "strict")
                            shortened_data = shorten_multiple_epoch(preprocessed_data, alignment_type, obj.event_duration);
                        elseif strcmp(obj.timing_mode, "dynamic")
                            shortened_data = preprocessed_data; % to be shortened later in edb generation
                        end
                        data_per_electrode = cellfun(@squeeze, num2cell(shortened_data, [3, 2]), 'UniformOutput', false);

                        database_row_temp.(table_col) = {data_per_electrode};
                    end
                    database = stack_tables(database, database_row_temp);
                end
                save(fullfile(obj.times_path, 'database.mat'), 'database', '-v7.3');
            else
                fprintf('Loading generic database... ');
                database = load(fullfile(obj.times_path, 'database.mat')).database;
                fprintf('Done.\n');
            end

            if ~isfile(fullfile(obj.times_path, 'electrodes_database.mat'))
                electrodes_database = obj.electrodes_table;
                electrodes_database.feat_set = cell([height(electrodes_database), 1]);
                electrodes_database.article_cluster_name = cell([height(electrodes_database), 1]);
                sub_nums_vector = electrodes_database.subject;
                surf_shapes = {'s', 'v', 'z', 'a', 'i'};
                [~, ~, ic] = unique(sub_nums_vector);
                surf_shapes_col = surf_shapes(ic)';
                electrodes_database.surf_shape = surf_shapes_col;

                for electrodes_row_idx = 1:height(electrodes_database)
                    row_sub = electrodes_database.subject(electrodes_row_idx);
                    row_align = electrodes_database.type(electrodes_row_idx);
                    row_idx_LocalProcessed = electrodes_database.idx_LocalProcessed(electrodes_row_idx);
                    row_electrode_id = electrodes_database.electrode(electrodes_row_idx);
                    row_cluster_label = electrodes_database.cluster_name{electrodes_row_idx};

                    if isempty(row_cluster_label)
                        electrodes_database.cluster_name{electrodes_row_idx} = 'NaN';
                        electrodes_database.article_cluster_name{electrodes_row_idx} = 'NaN';
                    end

                    if ~isnan(row_cluster_label) & ~strcmp(row_cluster_label, 'NaN')
                        dynamic_row_idx = ismember(obj.dynamic_window_timings.kuzdeba_cluster_label, row_cluster_label);
                    end

                    if any(dynamic_row_idx)
                        dynamic_timings_row = obj.dynamic_window_timings(dynamic_row_idx, :);
                        article_cluster_label = dynamic_timings_row.article_cluster_label{1};

                        electrodes_database.article_cluster_name{electrodes_row_idx} = article_cluster_label;

                        if row_align == 2
                            database_row = database.stim_data_per_electrode{database.sub_num == row_sub};
                            alignment = 'stim';
                        elseif row_align == 1
                            database_row = database.onset_data_per_electrode{database.sub_num == row_sub};
                            alignment = 'onset';
                        end

                        shortened_data = shorten_electrode_epoch(database_row{row_idx_LocalProcessed}, ...
                                                                dynamic_timings_row.start_ms_adjusted, ...
                                                                dynamic_timings_row.end_ms_adjusted, ...
                                                                alignment);

                        electrodes_database.feat_set(electrodes_row_idx) = {create_features(obj, row_electrode_id, shortened_data, dynamic_timings_row)};
                    else
                        electrodes_database.feat_set(electrodes_row_idx) = cell(1);
                    end
                end
                edb_row_num = table([1:height(electrodes_database)]', 'VariableNames', {'edb_row_num'});
                node_num = discretize([1:height(electrodes_database)]', obj.number_of_nodes);
                node_num_col = table(node_num, 'VariableNames', {'node_num'});
                word_elect_performance = table(zeros([height(electrodes_database), 1]), 'VariableNames', {'word_accuracy_change_wo_electrode'});
                cons_elect_performance = table(zeros([height(electrodes_database), 1]), 'VariableNames', {'cons_accuracy_change_wo_electrode'});
                vowel_elect_performance = table(zeros([height(electrodes_database), 1]), 'VariableNames', {'vowel_accuracy_change_wo_electrode'});
                electrodes_database = [edb_row_num, electrodes_database, node_num_col, word_elect_performance, cons_elect_performance, vowel_elect_performance];

                save(fullfile(obj.times_path, 'electrodes_database.mat'), 'electrodes_database', '-v7.3');
            else
                fprintf('Loading electrodes database... ');
                electrodes_database = load(fullfile(obj.times_path, 'electrodes_database.mat')).electrodes_database;
                fprintf('Done.\n');
            end
        end
    end
end

