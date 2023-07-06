sub_num = 357; 
brain_plot_class_label = 'word_name'; 
alignment = 'onset';

addpath /project/busplab/software/ecog/util/
addpath /project/busplab/software/spm12
addpath /project/busplab/software/conn
addpath /project/busplab/software/display	
addpath /project/busplab/software/display/surf

threshd_auc_data_with_mni = threshd_auc_data; 

for threshd_data_row_idx = 1:height(threshd_auc_data)
    stim_data = threshd_auc_data.stim_sig_data{threshd_data_row_idx}; 
    
    stim_mni_coords_col = zeros([height(stim_data), 3]); 
    for stim_data_row = 1:height(stim_data)
        stim_sub_num = stim_data.sub_num(stim_data_row); 
        stim_e_id = stim_data.stim_e_id(stim_data_row);
        stim_mni_coords_col(stim_data_row, :) = get_MNI_coords(stim_sub_num, stim_e_id, 'stim');        
    end
    
    stim_data_with_mni = addvars(stim_data, stim_mni_coords_col, 'NewVariableNames', 'stim_mni_coords', 'Before', 'mean_auc');
    threshd_auc_data_with_mni.stim_sig_data{threshd_data_row_idx} = stim_data_with_mni; 
    
    onset_data = threshd_auc_data.onset_sig_data{threshd_data_row_idx};
    
    onset_mni_coords_col = zeros([height(onset_data), 3]); 
    for onset_data_row = 1:height(onset_data)
        onset_sub_num = onset_data.sub_num(onset_data_row); 
        onset_e_id = onset_data.onset_e_id(onset_data_row);
        onset_mni_coords_col(onset_data_row, :) = get_MNI_coords(onset_sub_num, onset_e_id, 'onset');        
    end
    
    onset_data_with_mni = addvars(onset_data, onset_mni_coords_col, 'NewVariableNames', 'onset_mni_coords', 'Before', 'mean_auc');
    threshd_auc_data_with_mni.onset_sig_data{threshd_data_row_idx} = onset_data_with_mni;     
end

% plot_sig_electrodes_EZmode(threshd_auc_data_with_mni, sub_num, brain_plot_class_label, alignment);
