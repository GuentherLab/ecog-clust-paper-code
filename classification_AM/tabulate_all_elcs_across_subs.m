 %%%% use this script to compile lists of electrodes and MNI coordinates from all subjects
 % lists of electrodes will include those which were not responsive or included in Scott's analysis
 %
 % this is a necessary step for creating figure 1 of Meier et al. 2023, which illustrates the complete ecog grids and sEEG locations

set_paths()

unprocessed_elc_dir = [SOURCEDATA_DIR filesep 'unprocessed_electrodes_by_sub'];
sub_list = [357; 362; 369; 372; 376];
filename_coda = '_contact_locations.csv';
savefile_name = [SOURCEDATA_DIR filesep 'elcs_unprocessed'];
sourcedata_note = ...
    'data compiled by tabulate_all_elcs_across_subs.m from /ecog_clust/data/unpressed_electrodes_by_sub';

nsubs = length(sub_list);
elcs_unprocessed = table;
for isub = 1:nsubs
    thissub = sub_list(isub);
    sub_filename = [unprocessed_elc_dir filesep num2str(thissub) filename_coda];
    thissub_table = readtable(sub_filename);
    thissub_table.subject(:,1) = thissub;
    elcs_unprocessed = [elcs_unprocessed; thissub_table];
end
n_elcs = height(elcs_unprocessed);

elcs_unprocessed.mni_coord = [elcs_unprocessed.MNIX_mm_, elcs_unprocessed.MNIY, elcs_unprocessed.MNIZ];

is_depth = contains(elcs_unprocessed.Label,'depth','IgnoreCase',true); 
is_surface = contains(elcs_unprocessed.Label,'grid','IgnoreCase',true) ...
    | contains(elcs_unprocessed.Label,'strip','IgnoreCase',true) ...
    | contains(elcs_unprocessed.Label,'right temporal pole','IgnoreCase',true) ;
elcs_unprocessed.depth_type = cell(n_elcs,1);
elcs_unprocessed.depth_type(is_depth) = {'depth'};
elcs_unprocessed.depth_type(is_surface) = {'surface'};

elcs_unprocessed.right_hemi = nan(n_elcs,1);
is_right = contains(elcs_unprocessed.Label,'right','IgnoreCase',true);
is_left = contains(elcs_unprocessed.Label,'left','IgnoreCase',true);
elcs_unprocessed.right_hemi(is_right) = 1; 

elcs_unprocessed = removevars(elcs_unprocessed, {'MNIX_mm_','MNIY','MNIZ'});
elcs_unprocessed = movevars(elcs_unprocessed,...
    {'subject','Label','Contact','depth_type','mni_coord'},'Before',1);


save(savefile_name,'elcs_unprocessed',"sourcedata_note")