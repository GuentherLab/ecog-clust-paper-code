%%% plot all electrode locations for all subjects
% for fig 1

clear

setelc_ops.ecog_elc_shape = 's'; 
setelc_ops.ecog_elc_color = [0.7 0 0.7]; 

setelc_ops.seeg_elc_shape = 's';
setelc_ops.seeg_elc_color = [0.5 1 0]; 


set_paths()
if ~exist('electrodes','var')
    load(ALL_ELC_DATA_FILE)
end
setelc_ops.output_dir = [DATA_DIR filesep 'figs' filesep 'brainplots_all_electrodes'];

%%%% load all electrodes (includes those not analyzed by Scott Kuzdeba)
unprocessed_elc_file = [SOURCEDATA_DIR filesep 'elcs_unprocessed.mat'];
load(unprocessed_elc_file, 'elcs_unprocessed')
elcs_unprocessed.type(:,1) = nan; %%%% stim- vs. voice-aligned.... not important for this plot
elcs_unprocessed.global_clust_num(:,1) = nan; %%%% Kuzdeba cluster... not important for this plot

surf_ops.lh_tempfile = [DATA_DIR filesep 'lh_mnitemp_to_delete_', strrep(datestr(datetime),':','-'), '.txt']; % use new filename each time
surf_ops.rh_tempfile = [DATA_DIR filesep 'rh_mnitemp_to_delete_', strrep(datestr(datetime),':','-'), '.txt']; % use new filename each time
surf_ops.surfshowDir = 'C:/docs/code/matlab/display/surf';
surf_ops.clusterkey_file = [SOURCEDATA_DIR filesep 'clusterkey.mat'];

surf_ops.scale = 0.8; % electrode marker size
surf_ops.crop_borders_top_bot_left_right = [230,1100,240,1500]; % if not empty, crop saved image borders at this pixel values
% surf_ops.viewpoint = 'lateral';
surf_ops.add_legend = 0; 
surf_ops.colormap = [0 0 0];
surf_ops.force_close = 1; % if true, close the surf_show window after plotting/saving image




%% no elcs left
elc_temp = table; % empty table
surf_ops.selected_hemi = 'left';
surf_ops.viewpoint = 'left'; 
surf_ops.img_savename = [setelc_ops.output_dir filesep 'no_electrodes_brainplot_' surf_ops.selected_hemi];
plot_electrodes_on_brain_AM2(elc_temp,surf_ops);

%% no elcs right
elc_temp = table; % empty table
surf_ops.selected_hemi = 'right'; 
surf_ops.viewpoint = 'right'; 
surf_ops.img_savename = [setelc_ops.output_dir filesep 'no_electrodes_brainplot_' surf_ops.selected_hemi];
plot_electrodes_on_brain_AM2(elc_temp,surf_ops);

%% all elcs left
setelc_ops.sub = []; 
surf_ops.selected_hemi = 'left'; 
[elc_temp, surf_ops] = set_elcs(elcs_unprocessed,setelc_ops,surf_ops); 
surftab = plot_electrodes_on_brain_AM2(elc_temp,surf_ops);

%% all elcs right
setelc_ops.sub = []; 
surf_ops.selected_hemi = 'right'; 
[elc_temp, surf_ops] = set_elcs(elcs_unprocessed,setelc_ops,surf_ops); 
surftab = plot_electrodes_on_brain_AM2(elc_temp,surf_ops);

%% sub 357 left
setelc_ops.sub = 357; 
surf_ops.selected_hemi = 'left'; 
[elc_temp, surf_ops] = set_elcs(elcs_unprocessed,setelc_ops,surf_ops); 
surftab = plot_electrodes_on_brain_AM2(elc_temp,surf_ops);

%% sub 362 left
setelc_ops.sub = 362; 
surf_ops.selected_hemi = 'left'; 
[elc_temp, surf_ops] = set_elcs(elcs_unprocessed,setelc_ops,surf_ops); 
surftab = plot_electrodes_on_brain_AM2(elc_temp,surf_ops);

%% sub 372 left
setelc_ops.sub = 372; 
surf_ops.selected_hemi = 'left'; 
[elc_temp, surf_ops] = set_elcs(elcs_unprocessed,setelc_ops,surf_ops); 
surftab = plot_electrodes_on_brain_AM2(elc_temp,surf_ops);

%% sub 369 right
setelc_ops.sub = 369; 
surf_ops.selected_hemi = 'right'; 
[elc_temp, surf_ops] = set_elcs(elcs_unprocessed,setelc_ops,surf_ops); 
surftab = plot_electrodes_on_brain_AM2(elc_temp,surf_ops);

%% sub 362 right
setelc_ops.sub = 362; 
surf_ops.selected_hemi = 'right'; 
[elc_temp, surf_ops] = set_elcs(elcs_unprocessed,setelc_ops,surf_ops); 
surftab = plot_electrodes_on_brain_AM2(elc_temp,surf_ops);

%% sub 376 right
setelc_ops.sub = 376; 
surf_ops.selected_hemi = 'right'; 
[elc_temp, surf_ops] = set_elcs(elcs_unprocessed,setelc_ops,surf_ops); 
surftab = plot_electrodes_on_brain_AM2(elc_temp,surf_ops);


%%
function [elc_temp, surf_ops] = set_elcs(electrodes,setelc_ops,surf_ops)
    if isfield(setelc_ops,'sub') && ~isempty(setelc_ops.sub)  % use specific sub
        elc_temp = electrodes(electrodes.subject==setelc_ops.sub,:);
    else     % use all subs
        elc_temp = electrodes;
    end

    surf_ops.color_values_override = nan(height(elc_temp),3);

    elc_temp.surf_shape(elc_temp.depth_type == "surface",1) = {setelc_ops.ecog_elc_shape};
    surf_ops.color_values_override(find(elc_temp.depth_type == "surface"),1:3) = repmat(setelc_ops.ecog_elc_color,nnz(elc_temp.depth_type == "surface"),1);

    elc_temp.surf_shape(elc_temp.depth_type == "depth") = {setelc_ops.seeg_elc_shape};
    surf_ops.color_values_override(elc_temp.depth_type == "depth",1:3) = repmat(setelc_ops.seeg_elc_color, nnz(elc_temp.depth_type == "depth"), 1);


    surf_ops.shapes = char(elc_temp.surf_shape);
    surf_ops.viewpoint = surf_ops.selected_hemi;

    surf_ops.img_savename = [setelc_ops.output_dir filesep 'sub_' num2str(setelc_ops.sub) '_brainplot_all_elcs_' surf_ops.selected_hemi];
end