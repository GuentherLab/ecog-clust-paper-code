% if ~exist('params_struct', 'var')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%    File Path to Data    %%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     data_path = '/projectnb/busplab/Experiments/ECoG_fMRI_RS/Experiments/ECoG_Preprocessed_LJ/'; 
%     ml_path = fullfile(data_path, 'MLAlgoData'); 
%     mRMR_path = fullfile(data_path, 'MLResults/mRMR/');
%     compiled_data_path = fullfile(data_path, 'MLCompiledData');
%  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%    Define Parameters    %%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     event_duration = 2000;
%     window = 50;
%     stride = 25; 
% 
%     topN_feat = 150;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%   Constructing Parameter Struct   %%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%     times_folder = sprintf('e%s_w%s_s%s', string(event_duration), string(window), string(stride));
%      
%     params_struct = struct();
%     params_struct.data_path = data_path; 
%     params_struct.ml_path = ml_path;
%     params_struct.compiled_data_path = compiled_data_path;
%     params_struct.times_folder = times_folder;
%     params_struct.event_duration = event_duration;
%     params_struct.window = window;
%     params_struct.stride = stride;
% end

stim_all_nodes_all_iters = table();
onset_all_nodes_all_iters = table();

auc_results_path = fullfile(p.grouping_path, p.current_group_value, 'AUC_Results');
auc_res_file_list = dir(fullfile(auc_results_path, '*.mat'));
for node_idx = 1:numel(auc_res_file_list)
    load(fullfile(auc_results_path, auc_res_file_list(node_idx).name));
    
    stim_all_nodes_all_iters = [stim_all_nodes_all_iters; sub_stim_e_id_auc_all_labels_all_iter];
    onset_all_nodes_all_iters = [onset_all_nodes_all_iters; sub_onset_e_id_auc_all_labels_all_iter];
    
end
