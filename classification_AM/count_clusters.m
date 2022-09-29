% output numbers of electrodes in each cluster in each subject

clear
load('/usr2/postdoc/amsmeier/ECoG_Preprocessed_AM/electrodes.mat')

elc = electrodes;
subs = unique(elc.subject);
nsubs = length(subs); 

% voice-aligned analysis
clustnames = {'PtM-r', 'PtM-s', 'ME-sb', 'ME-sn', 'AP-r', 'AP-s' };
nclusts = length(clustnames); 
onset_sites_per_subject = table(NaN(nclusts, nsubs), 'VariableNames', {'n_357___362___369___372___376'}, 'RowNames', clustnames);
for isub = 1:nsubs
    thissub = subs(isub);
    matchrow = elc.type == 1 & elc.subject == thissub & ~isempty(elc.cluster_name);
    for iclust = 1:nclusts
        thisclust = clustnames{iclust};
        onset_sites_per_subject.n_357___362___369___372___376(iclust, isub) = nnz(matchrow & strcmp(elc.cluster_name, thisclust));
    end
end
onset_sites_per_subject.total_sites = sum(onset_sites_per_subject.n_357___362___369___372___376,2);

% stim-aligned analysis
clustnames = {'ESP-r', 'ESP-s', 'PtM-r', 'ME-sb', 'ME-sn', 'AP-s' };
nclusts = length(clustnames); 
stim_sites_per_subject = table(NaN(nclusts, nsubs), 'VariableNames', {'n_357___362___369___372___376'}, 'RowNames', clustnames);
for isub = 1:nsubs
    thissub = subs(isub);
    matchrow = elc.type == 2 & elc.subject == thissub & ~isempty(elc.cluster_name);
    for iclust = 1:nclusts
        thisclust = clustnames{iclust};
        stim_sites_per_subject.n_357___362___369___372___376(iclust, isub) = nnz(matchrow & strcmp(elc.cluster_name, thisclust));
    end
end
stim_sites_per_subject.total_sites = sum(stim_sites_per_subject.n_357___362___369___372___376,2);    

onset_sites_per_subject
stim_sites_per_subject

    