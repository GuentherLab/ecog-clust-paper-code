 %%% create dummy data for debugging areal_phonunit_preference.m

function [elc_out, subdata] = gen_dummy_data_for_bootstrap_test(genops)

n_phonunits = length(genops.phonunit_names); 

load(genops.elc_table_filename,'electrodes'); 
elc_out = electrodes(electrodes.type==1,:); %%%% keep only speech onset-aligned elcs; the stim-aligned elcs have type==2
n_elcs = height(electrodes); 
nsubs = height(genops.ntrials_by_sub);

 for isub = 1:nsubs
     thissub = genops.ntrials_by_sub.sub(isub);
     elcmatch = find(elc_out.subject == thissub); 
     nelcs_thissub = length(elcmatch); 
     ntrials_thissub = genops.ntrials_by_sub.ntrials(isub);
     for ielc_within_sub = 1:nelcs_thissub
         irow = elcmatch(ielc_within_sub);
         for iphonunit = 1:n_phonunits
             this_phon = genops.phonunit_names{iphonunit};
             %%% generate random accuracies for trials/phonunits in which this electrode was left out
             elc_out{irow,[this_phon, '_acc_without_elc']} = {round(rand(1,ntrials_thissub) + genops.acc_offset)};
         end
     end
 end

 %%%% generate random trialwise accuracies for each subject for when all electrodes are included
 subdata = genops.ntrials_by_sub;
for iphonunit = 1:n_phonunits
    this_phon = genops.phonunit_names{iphonunit};
    for isub = 1:nsubs
        ntrials_thissub = genops.ntrials_by_sub.ntrials(isub);
        %%% generate random accuracies for trials/phonunits in which all electrodes were included
        subdata{isub,[this_phon, '_acc_all_elcs']} = {round(rand(1,ntrials_thissub) + genops.acc_offset + genops.all_elcs_acc_boost)};
        subdata{isub,[this_phon '_acc_overall']} = mean(subdata{isub,[this_phon, '_acc_all_elcs']}{:}); 
    end
end

