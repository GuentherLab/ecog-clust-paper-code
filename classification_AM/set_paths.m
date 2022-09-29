%%%% set machine-specific paths for accessing data remotely
%
% scripts and functions for ECoG-Clust project
%
% updated 2022/7/20 by AM

compname  = getenv('COMPUTERNAME');

switch compname
    case 'MSI' % Andrew Meier laptop
        ROOT_DIR = 'B:'; % mapped drive
        % % % % % the following paths contain bad version of binocdf
        rmpath('C:\docs\code\matlab\packages\fieldtrip-20210616\external\stats')
        rmpath('C:\docs\code\matlab\packages\spm12\external\fieldtrip\external\stats\')
    otherwise
        ROOT_DIR = ''; %%% for accessing SCC directly via SSH
end

if isfolder(ROOT_DIR) % if SCC is currently available
    DATA_DIR = [ROOT_DIR '/projectnb/busplab/Experiments/ECoG_Preprocessed_AM']; % save results here
    SOFTWARE_DIR = [ROOT_DIR, '/project/busplab/software/ecog']; 
elseif ~isfolder(ROOT_DIR)% if SCC is not currently connected; save to local drive
    DATA_DIR = 'C:/Users/amsme/Downloads/clust'; % save results here
    SOFTWARE_DIR = 'C:/docs/code/matlab/packages'; 
end

% add software packages
addpath([SOFTWARE_DIR, '/ecog_clust'])
addpath([SOFTWARE_DIR, '/ecog_clust/classification_AM'])
addpath([SOFTWARE_DIR, '/ecog_clust/data'])
addpath([SOFTWARE_DIR, '/ecog_clust/util'])
addpath([SOFTWARE_DIR, '/spm12' ])
addpath([SOFTWARE_DIR, '/conn' ])
addpath([SOFTWARE_DIR, '/display' ])
addpath([SOFTWARE_DIR, '/display/surf' ])
    
clear compname