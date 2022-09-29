% gather results from stmf_toplevel, assign cluster labels
%
% updated by Andrew Meier 2021/5/25

vardefault('savedir', '/projectnb/busplab/Experiments/ECoG_Preprocessed_AM/stmf_results');
vardefault('prependstr', 'stmf_results_'); 

% vardefault('feature_to_analyze','word');
vardefault('feature_to_analyze','consonants');
% vardefault('feature_to_analyze','vowel');

vardefault('ops',struct);
field_default('ops', 'save_prepend', prependstr); 
filestr = [ops.save_prepend, feature_to_analyze]; 
dd = struct2table(dir(savedir));
stmfres_files = dd.name(contains(dd.name, filestr));
nfiles = size(stmfres_files,1); 

% concatenate data from all subjects
load([savedir, filesep, stmfres_files{1}], 'stmf_results', 'analysis_ops')
nsubs = length(analysis_ops.subjects_to_analyze); 
for ifile = 2:nfiles
    tmp = load([savedir, filesep, stmfres_files{ifile}], 'stmf_results'); 
    for isub = 1:nsubs
        thissub = analysis_ops.subjects_to_analyze(isub);
        % append results from this file to overall list
        stmf_results.trials{isub}.correct_ftrval = [stmf_results.trials{isub}.correct_ftrval, tmp.stmf_results.trials{isub}.correct_ftrval];
        stmf_results.iter_results{isub} = [stmf_results.iter_results{isub}; tmp.stmf_results.iter_results{isub}];
    end
end

% for each site, find all iterations in which it was included
% then give proportion correct for each site
%%% then create table which includes all sites across all subjects
sites = table; 
for isub = 1:nsubs
    isub
    nsites = height(stmf_results.sites{isub}); 
    niters = height(stmf_results.iter_results{isub});
    stmf_results.sites{isub}.iters_included = cell(nsites,1); % list of included iterations for each site
    stmf_results.sites{isub}.correct_prop = NaN(nsites,1); % avg proportion correct across all trials and included iterations
    for isite = 1:nsites
        site_idx = stmf_results.sites{isub}.idx_LocalProcessed(isite); % row within preprocessed_data.data
        this_site_iters = NaN(niters, 1); % preallocate for speed
        counter = 1; 
        for iter = 1:niters
            if ismember(site_idx, stmf_results.iter_results{isub}.idx_LocalProcessed{iter}) % if site was in this iter
                this_site_iters(counter) = iter; % add iter to this site's list
                counter = counter+1; 
            end
        end
        this_site_iters = this_site_iters(~isnan(this_site_iters)); 
        stmf_results.sites{isub}.iters_included{isite} = this_site_iters; 
        performlist = stmf_results.iter_results{isub}.frac_correct(this_site_iters); % performance across iters
        stmf_results.sites{isub}.correct_prop(isite) = mean(performlist); % average prop correct in all iters this site participated in
    end
    stmf_results.sites{isub}.zcorrect = zscore(stmf_results.sites{isub}.correct_prop); % zscore of site performance within subject
    stmf_results.sites{isub} = movevars(stmf_results.sites{isub}, {'zcorrect', 'correct_prop'}, 'After', 'cluster_name'); 
    sites = [sites; stmf_results.sites{isub}]; % add sites for this subject to the cross-subject table
end

parentdir = fileparts(savedir);
savestr = [parentdir, filesep, filestr];
save(savestr, 'stmf_results', 'sites', 'ops') % save compiled results
        
        