% clc

edb_var_info = who('-file', fullfile(p.times_path, 'electrodes_database.mat'));

if ismember('edb_cat', edb_var_info) && ismember('edb_nonan', edb_var_info)
    load(fullfile(p.times_path, 'electrodes_database.mat'), 'edb_cat', 'edb_nonan');
    
else
    p.current_group_value = 'none';
    files = dir(fullfile(p.grouping_path, p.current_group_value, 'excluded_electrode_data', '*.mat')); 

    filename = ['avg_acc_table_', p.times_label, '.mat']; 
    load(fullfile(p.grouping_path, p.current_group_value, filename));
    if ~exist('edb_mat', 'var')
        edb_mat = matfile(fullfile(p.times_path, 'electrodes_database.mat'), 'Writable', true); 
        edb_cat = edb_mat.electrodes_database; 
    end
        
    num_files = length(files); 

    for file_idx = 1:num_files
        file_name = files(file_idx).name;

        [row_num_start, row_num_end] = regexp(file_name, '(?<=row)\d+');
        row_num = str2double(file_name(row_num_start:row_num_end));

        sub_num = edb_cat.subject(row_num);

        [align_start, align_end] = regexp(file_name, '(stim|onset)');
        alignment_str = file_name(align_start:align_end);
        switch alignment_str
            case 'stim'
                alignment_type = 2;
            case 'onset'
                alignment_type = 1;
        end

        row_data = load(fullfile(files(file_idx).folder, files(file_idx).name));

        word_baseline_acc = avg_acc_table.word_name{avg_acc_table.sub_num == sub_num, alignment_type};
        cons_baseline_acc = avg_acc_table.consonants_name{avg_acc_table.sub_num == sub_num, alignment_type};
        vowel_baseline_acc = avg_acc_table.vowel_name{avg_acc_table.sub_num == sub_num, alignment_type};

        word_mean_acc = row_data.word_name;
        cons_mean_acc = row_data.consonants_name;
        vowel_mean_acc = row_data.vowel_name;

        edb_cat.word_accuracy_change_wo_electrode(row_num) = word_baseline_acc - word_mean_acc;
        edb_cat.cons_accuracy_change_wo_electrode(row_num) = cons_baseline_acc - cons_mean_acc;
        edb_cat.vowel_accuracy_change_wo_electrode(row_num) = vowel_baseline_acc - vowel_mean_acc;

    end
    

    % categorical clusters
    nan_cluster_order = {'ESP-s', 'ESP-r', 'PtM-s', 'PtM-r', 'ME-sb', 'ME-sn', 'AP-r', 'AP-s', 'NaN'};
    nonan_cluster_order = {'ESP-s', 'ESP-r', 'PtM-s', 'PtM-r', 'ME-sb', 'ME-sn', 'AP-r', 'AP-s'};

    edb_cat.cluster_name = categorical(edb_cat.cluster_name);
    edb_cat.cluster_name = reordercats(edb_cat.cluster_name, nan_cluster_order);
    
    % remove NaNs
    edb_nonan = edb_cat(~(edb_cat.cluster_name == 'NaN'), :);    
    edb_nonan.cluster_name = categorical(edb_nonan.cluster_name);
    edb_nonan.cluster_name = removecats(edb_nonan.cluster_name, 'NaN');
    edb_nonan.cluster_name = reordercats(edb_nonan.cluster_name, nonan_cluster_order);
    
    % save data to edb file 
    edb_mat.edb_cat = edb_cat;
    edb_mat.edb_nonan = edb_nonan; 
    
    save(fullfile(p.times_path, 'electrodes_database.mat'), 'edb_cat', 'edb_nonan', '-append');

end
