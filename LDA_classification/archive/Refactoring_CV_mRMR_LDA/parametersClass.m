classdef parametersClass 
    % parametersClass defines container of relevant info for analysis parameters
    %   Pass: (ecog_preprocessed_path, event_duration, window, stride,
    %   grouping_var, topN_feat_pool, topN_feat_indiv). 
    %   Outputs an object with params accessible through dot indexing. 

    properties(SetAccess = immutable)
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

        event_duration % The event duration specifies how much of the epoched signal will be used in the discretization process
        window % The window size (ms) of the moving average filter
        stride % The stride length (ms) of the MVA window

        topN_feat_pool % Top N features to select from the mRMR results in the pooled data analyses
        topN_feat_indiv % Top N features to select from the mRMR results in the individual data analyses
    end
    
    properties
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
    end
    
    methods 
        % Constructor:
        function params = parametersClass(ecog_preprocessed_path,...
                event_duration, window, stride, grouping_var,...
                topN_feat_pool, topN_feat_indiv)
            % Pass: (ecog_preprocessed_path, event_duration, window, stride,
            %     grouping_var, topN_feat_pool, topN_feat_indiv). 
            
            params.ecog_preprocessed_path = ecog_preprocessed_path; 
            params.epoched_data_path = fullfile(ecog_preprocessed_path, 'LocalEpoched');
            params.ml_path = fullfile(ecog_preprocessed_path, 'MachineLearning');
            params.algo_data_path = fullfile(params.ml_path, 'AlgoData');
            params.LDA_path = fullfile(params.ml_path, 'LDA');
            params.figs_quickref_path = fullfile(params.LDA_path, 'figs_quickref');
            params.times_label = sprintf('e%d_w%d_s%d', event_duration, window, stride);
            params.times_path = fullfile(params.LDA_path, params.times_label); 
            params.grouping_label = sprintf('grouped_by_%s', grouping_var);
            params.grouping_path = fullfile(params.times_path, params.grouping_label);
            params.full_run_label = strcat(params.times_label, '_', params.grouping_label);
            params.grouping_variable = grouping_var;
            params.event_duration = event_duration; 
            params.window = window; 
            params.stride = stride; 
            params.topN_feat_pool = topN_feat_pool; 
            params.topN_feat_indiv = topN_feat_indiv; 
            
            params.electrodes_table = load('electrodes.mat').electrodes; 
        end
        
        function [database, electrodes_database] = generate_database(obj)
            % generate_database generates database of subject level data, database of electrode level data.
            %    p = parametersClass(args) 
            %    [db, edb] = p.generate_database
            %    db is a table:
%
%                 sub_num    stim_data_per_electrode    onset_data_per_electrode
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

                    stim_electrodes_table = obj.electrodes_table(obj.electrodes_table.subject == sub_num & obj.electrodes_table.type == 2, :);
                    stim_electrodes_list = unique(stim_electrodes_table.idx_LocalProcessed);
%                     stim_electrodes_table = stim_electrodes_table(stim_electrodes_list, :);
                    
                    onset_electrodes_table = obj.electrodes_table(obj.electrodes_table.subject == sub_num & obj.electrodes_table.type == 1, :);
                    onset_electrodes_list = unique(onset_electrodes_table.idx_LocalProcessed);
%                     onset_electrodes_table = onset_electrodes_table(onset_electrodes_list, :);
                    
                    for file_idx = 1:numel(raw_data_files)
                        if contains(raw_data_filelist(file_idx).name, 'stim')
                            stim_AnalysisClass = load(fullfile(raw_data_filelist(file_idx).folder, raw_data_filelist(file_idx).name), 'preprocessed_data');
                            stim_preprocessed_data = stim_AnalysisClass.preprocessed_data.data; 
                            stim_shortened_data = shorten_multiple_epoch(stim_preprocessed_data, 'stim', obj.event_duration);
%                             stim_shortened_data = stim_shortened_data(stim_electrodes_list,:,:);
                            stim_data_per_electrode = cellfun(@squeeze, num2cell(stim_shortened_data, [3, 2]), 'UniformOutput', false); 
                        elseif contains(raw_data_filelist(file_idx).name, 'onset')
                            onset_AnalysisClass = load(fullfile(raw_data_filelist(file_idx).folder, raw_data_filelist(file_idx).name), 'preprocessed_data');
                            onset_preprocessed_data = onset_AnalysisClass.preprocessed_data.data;
                            onset_shortened_data = shorten_multiple_epoch(onset_preprocessed_data, 'onset', obj.event_duration);
%                             onset_shortened_data = onset_shortened_data(onset_electrodes_list,:,:);
                            onset_data_per_electrode = cellfun(@squeeze, num2cell(onset_shortened_data, [3, 2]), 'UniformOutput', false);
                        end
                    end

                    database_row_temp = table(sub_num, {stim_data_per_electrode}, {onset_data_per_electrode},...
                        'VariableNames', {'sub_num', 'stim_data_per_electrode', 'onset_data_per_electrode'});
                    database = [database; database_row_temp];

                    save(fullfile(obj.times_path, 'database.mat'), 'database', '-v7.3');
                end
            else 
                fprintf('Loading generic database... ');
                database = load(fullfile(obj.times_path, 'database.mat')).database;
                fprintf('Done.\n');                
            end
            
            if ~isfile(fullfile(obj.times_path, 'electrodes_database.mat'))
                
                electrodes_database = obj.electrodes_table; 
                electrodes_database.feat_set = cell([height(electrodes_database), 1]);
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
                    if row_align == 2
                        database_row = database.stim_data_per_electrode{database.sub_num == row_sub};
                    elseif row_align == 1
                        database_row = database.onset_data_per_electrode{database.sub_num == row_sub};
                    end

                    electrodes_database.feat_set(electrodes_row_idx) = {create_features(obj, row_electrode_id, database_row{row_idx_LocalProcessed})};

                    if isempty(electrodes_database.cluster_name{electrodes_row_idx})
                        electrodes_database.cluster_name{electrodes_row_idx} = 'NaN';
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

function feat_table = create_features(params, electrode_num, short_array)
    % Helper function for generate_database
sample_freq = 1000; 
window_samples = sample_freq * (params.window / 1000);
stride_samples = sample_freq * (params.stride / 1000);

if ndims(short_array) == 2
    movmean_array = movmean(short_array, window_samples, 2, 'Endpoints', 'discard');
    stride_idxs = 1:stride_samples:size(movmean_array, 2);
    movmean_array = movmean_array(:, stride_idxs); 
end

feat_labels_list = {};
for window_idx = 1:size(movmean_array, 2)
    feat_labels_list{window_idx} = sprintf('e%dw%d', electrode_num, window_idx);
end

feat_table = array2table(movmean_array, 'VariableNames', feat_labels_list);

end
