% defines parameter object containing relevant info for the pipeline

classdef parametersClass 
    properties(SetAccess = immutable)
        ecog_preprocessed_path
        epoched_data_path
        ml_path
        compiled_data_path
        mRMR_path
        times_path
        run_results_path
        times_flag
        full_run_label
        electrodes_table
        event_duration
        window
        stride
        topN_feat_pool
        topN_feat_indiv
        
    end
    properties (Constant)
        p_value = 0.05; 
        number_of_nodes = 50;        
        class_names = {'word_name', 'consonants_name', 'vowel_name'};
        class_names_formal = {'Word Name', 'Consonants Name', 'Vowel Name'};
        sub_nums = [357, 362, 369, 372, 376];
        sub_nums_formal = {'S357', 'S362', 'S369', 'S372', 'S376'};
        
    end
    
    methods 
        % Constructor:
        function params = parametersClass(ecog_preprocessed_path,...
                event_duration, window, stride, grouping_var,...
                topN_feat_pool, topN_feat_indiv)
                        
            params.ecog_preprocessed_path = ecog_preprocessed_path; 
            params.epoched_data_path = fullfile(ecog_preprocessed_path, 'LocalEpoched');
            params.ml_path = fullfile(ecog_preprocessed_path, 'MLAlgoData');
            params.compiled_data_path = fullfile(ecog_preprocessed_path, 'MLCompiledData');
            params.mRMR_path = fullfile(ecog_preprocessed_path, 'MLResults/mRMR');
            params.times_flag = sprintf('e%d_w%d_s%d', event_duration, window, stride);
            params.full_run_label = sprintf('e%d_w%d_s%d_grouped_by_%s', event_duration, window, stride, grouping_var);
            params.times_path = fullfile(params.mRMR_path, string(params.times_flag));
            params.run_results_path = fullfile(params.times_path, params.full_run_label);
            params.event_duration = event_duration; 
            params.window = window; 
            params.stride = stride; 
            params.topN_feat_pool = topN_feat_pool; 
            params.topN_feat_indiv = topN_feat_indiv; 
            
            load('electrodes.mat', 'electrodes');
            params.electrodes_table = electrodes; 
            
        end
    end
end


    
    
    

