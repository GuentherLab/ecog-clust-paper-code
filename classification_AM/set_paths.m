%%%% set machine-specific paths for accessing data remotely
%
% scripts and functions for ECoG-Clust project
%
% updated 2023/4/20 by AM

compname  = getenv('COMPUTERNAME');

switch compname
    case 'MSI' % Andrew Meier laptop
        ROOT_DIR = 'B:'; % mapped drive
        % % % % % the following paths contain bad version of binocdf
        rmpath('C:\docs\code\matlab\fieldtrip-20210616\external\stats')
        rmpath('C:\docs\code\matlab\spm12\external\fieldtrip\external\stats\')
    otherwise
        ROOT_DIR = ''; %%% for accessing SCC directly via SSH
end

if isfolder(ROOT_DIR) % if SCC is currently available
    DATA_DIR = [ROOT_DIR '/projectnb/busplab/Experiments/ECoG_Preprocessed_AM']; % save results here
    SOFTWARE_DIR = [ROOT_DIR, '/project/busplab/software/ecog']; 
elseif ~isfolder(ROOT_DIR)% if SCC is not currently connected; save to local drive
    DATA_DIR = 'C:/Users/amsme/Downloads/clust'; % save results here
    SOFTWARE_DIR = 'C:/docs/code/matlab'; 
end

SOURCEDATA_DIR = [SOFTWARE_DIR '\ecog_clust\data']; % load data from here
TOP_ELC_DATA_FILE = [SOURCEDATA_DIR filesep 'topelc_data_to_surf.mat'];
ALL_ELC_DATA_FILE = [SOURCEDATA_DIR '\electrodes.mat'];
FIG_DIR = [DATA_DIR, filesep, 'figs'];
DIR_LDA_2SEC_WINDOW = [DATA_DIR, filesep, 'LDA_outputs_2sec_window']; 
DIR_LDA_CLUSTER_SPEC_WINDOW = [DATA_DIR, filesep, 'LDA_outputs_cluster-specific_window']; 
TRIAL_ACC_FILENAME_2SEC_WINDOW = [DIR_LDA_2SEC_WINDOW filesep 'electrode_trialwise_accuracy.mat'];
TRIAL_ACC_FILENAME_CLUSTER_SPEC_WINDOW = [DIR_LDA_CLUSTER_SPEC_WINDOW filesep 'electrode_trialwise_accuracy.mat'];
SUBJECT_TABLE_FILENAME = [SOURCEDATA_DIR '\subjects_table.xlsx'];

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