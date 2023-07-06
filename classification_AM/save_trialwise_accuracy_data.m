%% current version of this script is WRONG; it uses data from when a single electrode was INCLUDED, not when it was EXCLUDED

%%%% reorganize accuracy from the leave-one-electrode-out analysis and save it
 % this script should be run before areal_phonunit_preference in order to organize necessary data

 clear
set_paths();

%%%%%%%%%%% the following option will determine which version of the machine learning analysis data to use
analysis_window_version = '2sec';
% analysis_window_version = 'cluster-specific';

subdata = table([357; 362; 369; 372; 376],[108; 143; 118; 144; 133],'VariableNames',{'sub','ntrials'});

path_string_elcdata = ['_name', filesep, 'sub_indiv_elect_cv_mRMR_LDA_data.mat']; 
path_string_subdata = ['_name', filesep, 'sub_cv_mRMR_LDA_data.mat']; 

% first line = shortened form to use in table variables
%%% second line = longer version used in Liam's saved data
phonunit_names = {'cons','vow','syl';'consonants','vowel','word'};

switch analysis_window_version
    case '2sec'
        LDA_data_dir = DIR_LDA_2SEC_WINDOW; 
        savefile = TRIAL_ACC_FILENAME_2SEC_WINDOW; 
    case 'cluster-specific'
        LDA_data_dir = DIR_LDA_CLUSTER_SPEC_WINDOW; 
        savefile = TRIAL_ACC_FILENAME_CLUSTER_SPEC_WINDOW;
end

% load electrode data; trialwise data will be filled into this table
load(ALL_ELC_DATA_FILE,'electrodes'); 
elc = electrodes(electrodes.type==1,:); %%%% keep only speech onset-aligned elcs; the stim-aligned elcs have type==2
clear electrodes
nelcs = height(elc); 


nsubs = height(subdata); 
for iphonunit = 1:size(phonunit_names,2)
    thisphon = phonunit_names{1,iphonunit}
    thisphon_long = phonunit_names{2,iphonunit};
    elc.([thisphon,'_acc_without_elc']) = cell(nelcs,1);
    subdata.([thisphon,'_acc_all_elcs']) = cell(nsubs,1);
    subdata.([thisphon,'_acc_overall']) = nan(nsubs,1);
    fname_to_load = [LDA_data_dir, filesep, thisphon_long, path_string_elcdata]; 
    load(fname_to_load)
    fname_to_load = [LDA_data_dir, filesep, thisphon_long,  path_string_subdata];
    load(fname_to_load)    
    
    for isub = 1:nsubs
        thissub = subdata.sub(isub);
        thissub_indiv_elcs = sub_indiv_elect_cv_mRMR_LDA_data.(['S',num2str(thissub)]).onset_data; 
        nfolds = height(thissub_indiv_elcs.onset_fold_label_table{1});
        nelcs_thissub = height(thissub_indiv_elcs);

        % individual electrode leave-one-out accuracy results
        for ielc = 1:nelcs_thissub
            elc_rowmatch = elc.subject==thissub & elc.electrode==thissub_indiv_elcs.onset_e_id(ielc); % this elc's row in all-subjects electrode table
            %%%% each row of the following table will be a data fold for when this electrode was left out
            elcdat = thissub_indiv_elcs.onset_fold_label_table{ielc}.([thisphon_long,'_name']);
            predicted_answers = vertcat(elcdat.onset_y_pred{:}); %%%% classifier prediction when this electrode was left out
            correct_answers = vertcat(elcdat.onset_y_test{:}); %%% actual value of the phonemic unit the classifier was trying to predict
            elc{elc_rowmatch,[thisphon,'_acc_without_elc']}{1} = [strcmp(predicted_answers,correct_answers)]'; 
        end

       % accuracy results when all elcs are included 
       thissub_all_elcs = sub_cv_mRMR_LDA_data.(['S',num2str(thissub)]).onset_data.([thisphon_long, '_name']);
       predicted_answers = vertcat(thissub_all_elcs.onset_y_pred{:});  %%%% classifier prediction using all els from this sub
       correct_answers = vertcat(thissub_all_elcs.onset_y_test{:}); %%% actual value of the phonemic unit the classifier was trying to predict
       subdata{isub,[thisphon '_acc_all_elcs']}  = {[strcmp(predicted_answers,correct_answers)]'}; 
       subdata{isub,[thisphon '_acc_overall']} = mean(subdata{isub,[thisphon, '_acc_all_elcs']}{:}); 
    end
end

save(savefile, 'elc', 'subdata', 'analysis_window_version')




 