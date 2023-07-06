%{

Final Pipeline for LDA Machine Learning Model with ECoG Data

Overview:
- Takes epoched ECoG data from LocalProcessed and constructs data structures
for following processing. 

- Creates feature sets according to user-defined temporal parameters. 

- For each alignment condition (visual stimulus and voice onset) and 
class label (word_name, consonants_name, vowel_name), feature sets are 
processed according to 4 paradigms: 
        All Electrodes (Categorical Labels)
        All Electrodes (One Hot Encoded Labels)
        Individual Electrodes (Categorical Labels)
        Individual Electrodes (One Hot Encoded Labels)

- In general, the processing for each paradigm follows this workflow:
    1. Split feature set into 5 cross validation folds. 
    2. Perform mRMR Feature Selection on 4 of 5 folds, denoted as the training
    set. 
    3. Use the topN features' data to train an LDA classifier. 
    4. Test classifier on remaining fold. 
    5. Repeat 2-4 using a different fold as the test set. 

- Compute Statistics on the results from above and present graphically. 

Contact Liam Jackson with questions: 
    lcj126@bu.edu
    413-579-4373

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%    Data Consolidation   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% This takes a while. Open the script to configure the pathing.  
MLData_Consolidation_final 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%    File Path to Data    %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_path = '/projectnb/busplab/Experiments/ECoG_fMRI_RS/Experiments/ECoG_Preprocessed_LJ/'; 
ml_path = fullfile(data_path, 'MLAlgoData'); 
mRMR_path = fullfile(data_path, 'MLResults/mRMR/');
compiled_data_path = fullfile(data_path, 'MLCompiledData');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%    Define Parameters    %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

class_labels = [{'word_name'}, {'consonants_name'}, {'vowel_name'}];

event_duration = 2000;
window = 50;
stride = 25; 

topN_feat = 150; % top 150 features from all data mRMR
topN_indiv_e = 30; % top 30 features per electrode mRMR
number_of_nodes = 50; % permutation test nodes
p_value = 0.05; % Z-scoring and Bonferroni correction

times_folder = sprintf('e%s_w%s_s%s', string(event_duration), string(window), string(stride));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%   Constructing Parameter Struct   %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params_struct = struct();
%%%% Pathes
params_struct.data_path = data_path; 
params_struct.ml_path = ml_path;
params_struct.mRMR_path = mRMR_path;
params_struct.times_folder = times_folder;

%%%% Data
params_struct.class_labels = class_labels;
params_struct.event_duration = event_duration;
params_struct.window = window;
params_struct.stride = stride;
params_struct.topN_feat = topN_feat;
params_struct.topN_indiv_e = topN_indiv_e;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  Create/Load Feature Set  %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%% Run this only if you haven't before:
fprintf('Creating Feature sets for data from all electrodes...\n');

create_mean_HG_features_reg_ohe_standalone;

fprintf('Done. Now creating Feature sets for data from individual electrodes...\n');

indiv_electrode_data_preparation_reg_ohe;

fprintf('Done.\n');

%%%% Otherwise just load the feature set here:
fprintf('Loading feature sets for all electrodes combined... ');

feat_set_folder = fullfile(compiled_data_path, 'FeatureSets', times_folder);
load(fullfile(feat_set_folder, 'feature_set_table'), 'all_subs_feat_table');

fprintf('Done. \nNow loading feature sets for individual electrodes... ');

indiv_elect_feat_set_folder = fullfile(compiled_data_path, 'FeatureSets', times_folder);
load(fullfile(indiv_elect_feat_set_folder, 'indiv_elect_feature_set_table'), 'all_subs_indiv_elect_feat_table');

fprintf('Done.\n');

params_struct.all_subs_feat_table = all_subs_feat_table;
params_struct.all_subs_indiv_elect_feat_table = all_subs_indiv_elect_feat_table;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   Begin the Training/Analysis   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
all_class_sub_cv_mRMR_LDA_OvR_data = struct();
all_class_sub_cv_mRMR_LDA_data = struct();
all_class_sub_indiv_elect_cv_mRMR_LDA_OvR_data = struct();
all_class_sub_indiv_elect_cv_mRMR_LDA_data = struct();

class_labels = [{'word_name'}, {'consonants_name'}, {'vowel_name'}];
for class_label_idx = 1:length(class_labels)
    class_label = class_labels{class_label_idx};
    params_struct.class_label = class_label;
    fprintf('Class label: %s\n', class_label);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%% Step 1. all_electrodes_cv_mRMR_LDA_OvR.m
    fprintf('Step 1. All Electrodes, One vs. Rest\n');
%%%%%%%% Only run this if you have to:
    sub_cv_mRMR_LDA_OvR_data = all_electrodes_cv_mRMR_LDA_OvR(params_struct);
%_________________________________________________________________________%
%%%%%%%% Otherwise, run this to load the data you already have:
%     indiv_elect_cv_mRMR_LDA_OvR_res_file = sprintf('e%s_w%s_s%s_sub_indiv_elect_cv_mRMR_LDA_OvR_data', string(event_duration), string(window), string(stride));
%     indiv_elect_cv_mRMR_LDA_OvR_res_filename_full = fullfile(mRMR_path, times_folder, class_label, indiv_elect_cv_mRMR_LDA_OvR_res_file);
%     load(indiv_elect_cv_mRMR_LDA_OvR_res_filename_full, 'sub_indiv_elect_cv_mRMR_LDA_OvR_data');  

%_________________________________________________________________________%
    all_class_sub_cv_mRMR_LDA_OvR_data.(class_label) = sub_cv_mRMR_LDA_OvR_data;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 2. all_elect_cv_mRMR_LDA.m
    fprintf('Step 2. All Electrodes, Regular Labels\n');
%%%%%%%% Only run this if you have to:
    sub_cv_mRMR_LDA_data = all_electrodes_cv_mRMR_LDA(params_struct); 
%_________________________________________________________________________%    
%%%%%%%% Otherwise, run this to load the data you already have:
%     cv_mRMR_LDA_res_file = sprintf('e%s_w%s_s%s_sub_cv_mRMR_LDA_data', string(event_duration), string(window), string(stride));
%     cv_mRMR_LDA_res_filename_full = fullfile(mRMR_path, times_folder, class_label, cv_mRMR_LDA_res_file);
%     load(cv_mRMR_LDA_res_filename_full, 'sub_cv_mRMR_LDA_data'); 

%_________________________________________________________________________%    
    all_class_sub_cv_mRMR_LDA_data.(class_label) = sub_cv_mRMR_LDA_data;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 3. indiv_elect_cv_mRMR_LDA_OvR.m
    fprintf('Step 3. Individual Electrodes, One vs. Rest\n');
%%%%%%%% Only run this if you have to:
    sub_indiv_elect_cv_mRMR_LDA_OvR_data = indiv_electrode_cv_mRMR_LDA_OvR(params_struct);
%_________________________________________________________________________%
%%%%%%%% Otherwise, run this to load the data you already have:
%     indiv_elect_cv_mRMR_LDA_OvR_res_file = sprintf('e%s_w%s_s%s_sub_indiv_elect_cv_mRMR_LDA_OvR_data', string(event_duration), string(window), string(stride));
%     indiv_elect_cv_mRMR_LDA_OvR_res_filename_full = fullfile(mRMR_path, times_folder, class_label, indiv_elect_cv_mRMR_LDA_OvR_res_file);
%     load(indiv_elect_cv_mRMR_LDA_OvR_res_filename_full, 'sub_indiv_elect_cv_mRMR_LDA_OvR_data');  

%_________________________________________________________________________%    
    all_class_sub_indiv_elect_cv_mRMR_LDA_OvR_data.(class_label) = sub_indiv_elect_cv_mRMR_LDA_OvR_data;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 4. indiv_electrode_cv_mRMR_LDA.m
    fprintf('Step 4. Individual Electrodes, Regular Labels\n');
%%%%%%%% Only run this if you have to:
    sub_indiv_elect_cv_mRMR_LDA_data = indiv_electrode_cv_mRMR_LDA(params_struct);
%_________________________________________________________________________%
%%%%%%%% Otherwise, run this to load the data you already have:
%     indiv_elect_cv_mRMR_LDA_res_file = sprintf('e%s_w%s_s%s_sub_indiv_elect_cv_mRMR_LDA_data', string(event_duration), string(window), string(stride));
%     indiv_elect_cv_mRMR_LDA_res_filename_full = fullfile(mRMR_path, times_folder, class_label, indiv_elect_cv_mRMR_LDA_res_file);
%     load(indiv_elect_cv_mRMR_LDA_res_filename_full, 'sub_indiv_elect_cv_mRMR_LDA_data');  

%_________________________________________________________________________%        
    all_class_sub_indiv_elect_cv_mRMR_LDA_data.(class_label) = sub_indiv_elect_cv_mRMR_LDA_data; 
    
end

%%%% This is a useful block to have, next time you need it it won't take 15
%%%% years to run. Just load final_checkpoint.mat

checkpoint_file = sprintf('checkpoint_e%s_w%s_s%s.mat', string(event_duration), string(window), string(stride));
checkpoint_filename_full = fullfile(mRMR_path, times_folder, checkpoint_file);
save(checkpoint_filename_full,...
    'all_class_sub_cv_mRMR_LDA_OvR_data',...
    'all_class_sub_cv_mRMR_LDA_data',...
    'all_class_sub_indiv_elect_cv_mRMR_LDA_OvR_data',...
    'all_class_sub_indiv_elect_cv_mRMR_LDA_data',... 
    '-v7.3');
    
%%

% load(checkpoint_filename_full,...
%     'all_class_sub_cv_mRMR_LDA_OvR_data',...
%     'all_class_sub_cv_mRMR_LDA_data',...
%     'all_class_sub_indiv_elect_cv_mRMR_LDA_OvR_data',...
%     'all_class_sub_indiv_elect_cv_mRMR_LDA_data');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   Compiling Results   %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 5. indiv_electrode_ovr_auc_compiled.m
fprintf('Step 5. Calculating mean AUC scores and compiling data... ');
indiv_electrode_ovr_auc_compiled; 
fprintf('Done.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 6. (multinodal scripts)
%%%%%%%% Only run this if you have to:
%%%% WARNING: the following function, depending on the number of nodes, may
%%%%%%%% require multiple hours/days to complete.
% fprintf('Step 6. Now executing multinodal processing to create random data... ');
% job = execute_multinode_pipeline(number_of_nodes);
% fprintf('Done.\nNow compiling results... ');
% compile_all_nodes_all_iterations;
% fprintf('Done.\n\n');

%%%%%%%% Otherwise, run this to load the data you already have:

fprintf('Step 6. Now compiling multinodal results... ');
compile_all_nodes_all_iterations;
fprintf('Done.\n\n');

%%%% Thresholding real data with respect to significance indicated by
%%%% permutation test results
threshold_table = find_significant_auc_threshold(stim_all_nodes_all_iters, onset_all_nodes_all_iters, p_value);
threshd_auc_data = threshold_indiv_data(sub_stim_e_id_auc_all_labels, sub_onset_e_id_auc_all_labels, threshold_table);

close all; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   Results and Stats   %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figures_path = fullfile(mRMR_path, times_folder, 'Figures');
mkdir(figures_path);

%%%% Step 7. Visualizing results for all electrodes, regular labels
all_reg_visuals;        % fig 1
%%%% Step 8. Visualizing results for all electrodes, ovr labels
all_ovr_visuals;        % fig 2-4
%%%% Step 9. Visualizing distribution of aucs based on random data
perm_test_norm_visual;  % fig 5
%%%% Step 10. Visualizing step 8, but significant electrodes only
sig_ovr_visuals;        % fig 6-8
%%%% Step 11. Visualizing proportions of electrodes considered significant
ind_ovr_e_proportion;   % fig 9-11
%%%% Step 12. Projecting significant electrodes onto brain
mni_electrode_visual;   % fig 12 (not a figure, creates a surf_show plot)

close all; 




