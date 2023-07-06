 %%% determine which trials have analyzable ecog responses, then find the behavioral data for these trials
%
% called by stmf_classification
%
%% updated 2021/4/17 by Andrew Meier
function filelist = match_trial_to_behavior()

filelist_xlsx = 'filelist_for_classification.xlsx'; % list of relevant files for each subject
% ecog_localprocessed_topdir = fullfile(ecog_preprocessed_path, 'LocalProcessed'); % directory containing analyzed ecog trial data folders for each subject
target_words_file = 'Target_Words.mat'; % list of target words to be pronounced
stim_index_key_filename = 'stim_index_key.xlsx'; % numerical indices for stim features

nblocks = 4; % maximum number of block of experimental data per subject; must correspond to number of blocks in the xlsx filelist
max_trials_per_block = 36; % number of syllable pairs; we are only analyzing the first of each pair
vars_to_expand = {'analyze_this_block', 'ntrials_by_block', 'missing_trials_by_block', 'block_name'}; % variables to be formatted

%%% import and organize the file table
filelist = readtable(filelist_xlsx);
nsubs = height(filelist); 
[~,~,var_row] = xlsread(filelist_xlsx);     var_row = var_row(1,:); 
nameinds = cellfun(@isstr,var_row);
filelist.Properties.VariableNames(nameinds) =  var_row(nameinds);   %%% plug in variable names
% convert all strings to numbers
cell_var_inds = false(1,length(nameinds)); 
for icol = 1:size(filelist,2)
    if iscell(filelist{1,icol})  && ischar(filelist{1,icol}{1})
        cell_var_inds(icol) = true; 
        filelist{:,icol} = cellfun(@str2num,filelist{:,icol},'UniformOutput',false);
    end
end
temptable = filelist; 
for ivar = 1:length(vars_to_expand) % group all blocks under one table variable
    thisvar = vars_to_expand{ivar};
    old_table_ind = find(strcmp(thisvar,var_row));
    new_table_ind = find(strcmp(thisvar, temptable.Properties.VariableNames)); 
    cellvars_this_group = cell_var_inds(old_table_ind:old_table_ind+nblocks-1);
    if any(cellvars_this_group) %%% if any vars to be grouped are cells, make all cells
        varinds_to_make_cells = new_table_ind - 1 + find(~cellvars_this_group); 
        for i = 1:length(varinds_to_make_cells)
            newcol =  num2cell(temptable{:,varinds_to_make_cells(i)});
            eval(['temptable.', temptable.Properties.VariableNames{varinds_to_make_cells(i)}, '=newcol;'])
        end
    end
    temptable = mergevars(temptable, new_table_ind:[new_table_ind+nblocks-1], 'NewVariableName',thisvar);
end
filelist = temptable; clear temptable

%%% make table of target words
load(target_words_file);
cellcol = cell(max_trials_per_block,1); nancol = nan(max_trials_per_block,1); 
stim = table(Word(1:2:end), cellcol,                    cellcol,       nancol,    nancol,       nancol, 'VariableNames',...
            {   'word_name', 'consonants_name', 'vowel_name',  'word', 'consonants', 'vowel'}); %%% keep only the first produced word of each pair
chars = char(stim.word_name);
stim.consonants_name = string(chars(:,[1, 4])); % get the phonemes of each word
stim.vowel_name = string(chars(:,2));

%%% assign each stim feature an index
stim_index_key = readtable(stim_index_key_filename); % load preset value/index list
for istim = 1:max_trials_per_block
    stim.word(istim) = stim_index_key.index(strcmp(stim.word_name{istim}, stim_index_key.feature_val));
    stim.consonants(istim) = stim_index_key.index(strcmp(stim.consonants_name{istim}, stim_index_key.feature_val));
    stim.vowel(istim) = stim_index_key.index(strcmp(stim.vowel_name{istim}, stim_index_key.feature_val));
end

%%% get trials to use by excluding missing trials
%   usable trials are any that are not marked as NaN in subject's OnsetTable.mat within /projectnb/busplab/Experiments/ECoG_Preprocessed/LocalProcessed
%%%    the indices listed in usable_trials_by_block point to trial indices in the stim table
filelist.usable_trials_by_block = cell(size(filelist.missing_trials_by_block)); 
for iblock = 1:length(filelist.usable_trials_by_block(:))
    filelist.usable_trials_by_block{iblock} = 1:max_trials_per_block;
    if ~isempty(filelist.missing_trials_by_block{iblock}) &&  ~isnan(filelist.missing_trials_by_block{iblock}(1))  % if there are missing trials in this block
        missing_trials_this_block = filelist.missing_trials_by_block{iblock};
        filelist.usable_trials_by_block{iblock}(missing_trials_this_block) = [];  %%% exclude the missing trials from list of analyzable trials
    end
end

% match each usable trial to its stim features; find blocks to skip
filelist.trials = cell(nsubs, 1); %%% for each subject, create a table containing data about each usable trial
filelist.blocks_to_skip = cell(nsubs, 1);  %%% list of blocks to not analyze any trials from in subsequent analyses
for isub = 1:nsubs
    itrial_within_block = []; % initialize
    block = []; %%% initialize; this variable will label the block of each trial
    analyze = false(0,0); % initialize; this variable will indicate whether a trial should be analyzed (based on inclusion in good vs. bad blocks)
    for iblock = 1:nblocks
        itrial_within_block = [itrial_within_block; filelist.usable_trials_by_block{isub, iblock}']; % add the trial inds from this block, this subject
        block = [block; repmat(iblock, max([filelist.ntrials_by_block(isub,iblock), 0]), 1)]; %%% block labels to be assigned each trial within trials table
        analyze = [analyze; repmat(  logical(filelist.analyze_this_block(isub,iblock)), max([filelist.ntrials_by_block(isub,iblock), 0]), 1  )]; %  analyze these trials or not
    end
    trials = [table(block, itrial_within_block, analyze), stim([itrial_within_block], :)]; %%% get stim data for all trials for this sub
    filelist.trials{isub} = trials; clear itrial_within_block block analyze trials      %%% add trials data to filelist
end
clear nsubs
save('filelist.mat', 'filelist');

end
