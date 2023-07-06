function database = create_database(params)

sub_nums = params.sub_nums;
sub_nums_formal = params.sub_nums_formal;
database = table();
for sub_num_formal_idx = 1:length(sub_nums_formal)
    sub_num = sub_nums(sub_num_formal_idx);
    sub_num_formal = sub_nums_formal{sub_num_formal_idx};
    
    raw_data_filelist = dir([fullfile(params.epoched_data_path, sub_num_formal) '*/*.mat']);
    raw_data_files = {raw_data_filelist.name};
    
    for file_idx = 1:numel(raw_data_files)
        if contains(raw_data_filelist(file_idx).name, 'stim')
            stim_AnalysisClass = load(fullfile(raw_data_filelist(file_idx).folder, raw_data_filelist(file_idx).name), 'preprocessed_data');
            stim_preprocessed_data = stim_AnalysisClass.preprocessed_data.data; 
            stim_shortened_data = shorten_multiple_epoch(stim_preprocessed_data, 'stim', params.event_duration);
        elseif contains(raw_data_filelist(file_idx).name, 'onset')
            onset_AnalysisClass = load(fullfile(raw_data_filelist(file_idx).folder, raw_data_filelist(file_idx).name), 'preprocessed_data');
            onset_preprocessed_data = onset_AnalysisClass.preprocessed_data.data;
            onset_shortened_data = shorten_multiple_epoch(onset_preprocessed_data, 'onset', params.event_duration);
        end
    end
    
    database_row_temp = table({sub_num_formal}, sub_num,...
        table({stim_shortened_data}, {onset_shortened_data}, 'VariableNames', {'stim', 'onset'}),...
        'VariableNames', {'sub_num_formal', 'sub_num', 'raw_data'});
    database = [database; database_row_temp];

end





