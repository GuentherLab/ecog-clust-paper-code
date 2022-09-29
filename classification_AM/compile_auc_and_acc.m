 % compile the accuracy and area-under-the-curve (AUC) data created by Liam Jackson 
 % ...... match each electrode in the electrodes.mat table with its accuracy and AUC results
%
%%% 2022-2-2 by Andrew Meier

electrodes_file = '/projectnb2/busplab/Experiments/ECoG_Preprocessed_AM/electrodes.mat';
% the following contains folders for acc and auc results for word, consonant, vowel
acc_auc_parent_dir = '/projectnb2/busplab/Experiments/ECoG_Preprocessed_LJ/MachineLearning/LDA/e2000_w50_s25/grouped_by_none/none';
    acc_data_filename = 'sub_indiv_elect_cv_mRMR_LDA_data.mat'; 
    auc_data_filename = 'sub_indiv_elect_cv_mRMR_LDA_OvR_data.mat'; 
stim_key_file = '/projectnb2/busplab/Experiments/ECoG_Preprocessed_AM/stim_index_key.xlsx';
label_types = {'consonants', 'vowel', 'word'};
align_conds = {'onset', 'stim'}; 
sub_nums = [357; 362; 369; 372; 376]; 
%
savefile = '/usr2/postdoc/amsmeier/electrodes_acc_auc';



load(electrodes_file, 'electrodes')
nsubs = length(sub_nums);
sub_names = {};
for isub = 1:nsubs
    sub_names{isub} = ['S', num2str(sub_nums(isub))];
end
nancol = NaN(height(electrodes), 1); 
addtable = table(nancol, nancol, nancol, nancol, nancol, nancol, ...
    'VariableNames', {'acc_consonants', 'acc_vowel', 'acc_word', 'auc_consonants', 'auc_vowel', 'auc_word'}); 
electrodes = [electrodes, addtable]; 
stimkey = readtable(stim_key_file);


n_align_conds = length(align_conds); 
n_label_types = length(label_types);
for i_label_type = 1:n_label_types
    clear thissub this_align sub_indiv_elect_cv_mRMR_LDA_data sub_indiv_elect_cv_mRMR_LDA_OvR_data % clear previous data
    this_label_type = label_types{i_label_type}
    labelvals = stimkey.feature_val(strcmp(stimkey.type, this_label_type)); % list of vals for this label type
    n_labelvals = length(labelvals); 
    this_dir = [acc_auc_parent_dir, filesep, this_label_type, '_name'];

   
    load([this_dir, filesep, acc_data_filename]) % load accuracy data for this label type
    load([this_dir, filesep, auc_data_filename]) % load AUC data for this label type
    %%
    for isub = 1:nsubs
        this_sub_num = sub_nums(isub); 
        thissub = sub_names{isub};
        
        for i_align = 1:n_align_conds
            this_align = align_conds{i_align};
             % get number of folds that the data was divided into for machine learning cross-validation
            nfolds = height(sub_indiv_elect_cv_mRMR_LDA_data.(sub_names{1}).([this_align, '_data']).([this_align, '_fold_label_table']){1});
            if strcmp(this_align, 'onset')
                align_number = 1; % for comparing to electrodes table
            elseif strcmp(this_align, 'stim')
                align_number = 2; % for comparing to electrodes table
            else 
                error('bad alignment type')
            end
                
            % get number of electrodes in this subject and this align condition
            n_elc_this_sub_and_align = height(sub_indiv_elect_cv_mRMR_LDA_data.(thissub).([this_align, '_data']).([this_align, '_fold_label_table']));
            for ielc = 1:n_elc_this_sub_and_align
                % find corresponding electrodes table row
                this_electrode_ind = sub_indiv_elect_cv_mRMR_LDA_data.(thissub).([this_align, '_data']).([this_align, '_e_id'])(ielc);
                electrode_row = find([electrodes.subject == this_sub_num] & [electrodes.type == align_number] & ...
                    [electrodes.electrode == this_electrode_ind]);
                
                % get accuracy
                acc_all_folds = sub_indiv_elect_cv_mRMR_LDA_data.(thissub).([this_align, '_data']). ... % accuracy from all cross-validation folds
                      ([this_align, '_fold_label_table']){ielc}.([this_label_type, '_name']).([this_align, '_acc']);
                mean_acc = mean(acc_all_folds); 
                electrodes.(['acc_', this_label_type])(electrode_row) = mean_acc; 
        
                % get AUC
                clear auc_table
                auc_table = table(labelvals, NaN(n_labelvals, nfolds), NaN(n_labelvals, 1) , 'VariableNames', {'label', 'auc_all_folds', 'mean_auc'});
                for ival = 1:n_labelvals
                    thisval = labelvals{ival};
                    auc_table.auc_all_folds(ival,:) = sub_indiv_elect_cv_mRMR_LDA_OvR_data.(thissub).([this_align, '_data']). ... % accuracy from all cross-validation folds
                        ([this_align, '_fold_label_table']){ielc}.(thisval).([this_align, '_auc']);
                end
                auc_table.mean_auc = mean(auc_table.auc_all_folds,2); % get auc from all folds for each label
                mean_auc_all_labels = mean(auc_table.mean_auc); % mean auc across all labels for this electrode/align_condition/label_type
                electrodes.(['auc_', this_label_type])(electrode_row) = mean_auc_all_labels; 
            end
        end
        save(savefile,'electrodes') % save compiled results after each subject
    end
end
                
 
            
            
            
            
            
            