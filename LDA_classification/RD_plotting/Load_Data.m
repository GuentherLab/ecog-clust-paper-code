% sub_str_all = {'S357', 'S362', 'S369', 'S372', 'S376'};
% base_data_path = '/projectnb/busplab/Experiments/ECoG_Preprocessed_RD/LocalEpoched/'; 
% data_type = 'onset_12'; % 'onset_12', 'stimuli_12' etc. (file contains this string in name)
% 
% for sub_ind = 1:1%length(all_sub_str)
%     sub_str = sub_str_all{sub_ind};
%     
%     % Stuff
%     
% end

%%
function [event_info, preprocessed_data] = Load_Data(base_data_path, data_type, sub_str)
% Returns file containing data_type from path specified by base_data_path

    data_path = fullfile(base_data_path, sub_str);
    folder_contents = dir(data_path);
    
    % Find indices of folder that contain data_type in its name
    ids = zeros(1:length(folder_contents));
    for i = 3:length(folder_contents) % ignore '.' and '..'
       file_name = folder_contents(i).name;
       
       ids(i) = i*contains(file_name, data_type);
    end
    ids = ids(ids~=0);
    
    % Return whether or not a file was found
    if length(ids) ~= 1
        fprintf(['Verify base_data_path and data_type\n'...
            '%d files found in location\n'...
            '\s'], length(ids),fullfile(data_path));
    else
        load(fullfile(data_path, folder_contents(ids).name),...
            'event_info', 'preprocessed_data')
    end
end % End Load_Data