%%% this script gets all channel indices for each subject from LocalEpoched
%
% 'folder_full_path' = folder where preprocessed data will be loaded from... 
% ... probably /projectnb/busplab/Experiments/ECoG_Preprocessed/LocalEpoched
% 'align_condition' should be either 'onset' or 'stimuli'
%
% modified by Andrew Meier 2021/5/23 from Rohan's Generate_Chan_IDs_mat.m

function Generate_Chan_IDs_AM(folder_full_path, align_condition)

save_dir = '/projectnb/busplab/Experiments/ECoG_Preprocessed_AM/old'; % save output here

% No folder called
if ~nargin
    fprintf('Error: No file selected\n')
    return
end

% Check if you want to save
save_flag = true;

if save_flag
    temp = split(folder_full_path,'/');
    temp = temp(~cellfun('isempty',temp));
    temp{end} = strcat(temp{end}, '_', align_condition, '_chan_ids.mat');
    save_file_path = fullfile(save_dir, temp{end});
    if exist(save_file_path,'file')
        temp = input(['File exists at\n', save_file_path,...
            '\nContinue and overwrite? (y/n)?\n'],'s');
        if ~strcmpi(temp,'y')
            return
        end
    else
        fprintf(['Creating file at\n', save_file_path,'\n'])
    end
end

% Get patient folders
folder_info = dir(folder_full_path);
subfolder_flags = [folder_info.isdir];
subfolders = folder_info(subfolder_flags);

% Exclude '.' and '..' subfolders
temp = [];
for i = 1:numel(subfolders)
   if any(strcmp(subfolders(i).name, {'.', '..'}))
      temp(end+1) = i;
   end    
end
subfolders(temp) = [];

data_filename = ['Epoch_', align_condition, '_1_Hilbert_HG.mat']; % all data files to load will have this name

% Get chan_list from each subfolder and save to chan_list struct
chan_list = struct;
for i = 1:numel(subfolders)
   % Load .mat file and find preprocessed_data
   preproc_file = dir(fullfile(folder_full_path,subfolders(i).name,data_filename));
   load(fullfile(folder_full_path,subfolders(i).name,preproc_file(1).name),'preprocessed_data')
   
   % If preprocessed_data was in the .mat loaded, add its chan_list
   if ~exist('preprocessed_data','var')
       fprintf(['Error: preprocessed_data not found in:\n' ...
           fullfile(folder_full_path,subfolders(i).name,preproc_file(1).name)])
   else
       % Add preprocessed_data.chan_ids to return_variable
       chan_list.(subfolders(i).name) = preprocessed_data.chan_ids;
   end    
end

% Save
if save_flag
    save(save_file_path,'chan_list')    
end
end