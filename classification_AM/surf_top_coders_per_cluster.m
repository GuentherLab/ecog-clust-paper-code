 %%% plot 1 panel per cluster containing best encoders
 % top encoders will be drawn from cluster irrespective of cluster label.....
 % ...... so different clusters will have different proportions of top encoders
 %
 % this script only saves individual brain images; to put them into 1 figure, use compile_top_coders_surf.m
 % 
 %%% updated 2022/7/13 by AM

set(0,'DefaultFigureWindowStyle','normal')
%     set(0,'DefaultFigureWindowStyle','docked')

set_paths()

% data processed by the script load_top_encoders.m
topelc_data_filename = [ROOT_DIR, filesep, 'projectnb/busplab/Experiments/ECoG_Preprocessed_AM/topelc_data_to_surf'];

% plotting options
ops.subjects_as_shapes = 0;
ops.force_close = 1; 
ops.add_legend = 0; 
ops.add_title = 0; 
ops.clusters_as_colors = 0; 

cons_color = [1 0 0]; 
vow_color = [0 1 0];
syl_color = [0 0 1]; 

clusts_to_plot = 1:6; 
%     clusts_to_plot = 4;

savedir = [ROOT_DIR, filesep, 'projectnb/busplab/Experiments/ECoG_Preprocessed_AM/top_elcs_single_clust']; 
    
   

%%
load(topelc_data_filename)

for iclust = clusts_to_plot % 1:nclusts
    thisclust_global = global_clust_num_list(iclust); 
%     
% %     left hemisphere
   ops.selected_hemi = 'left'; 
   ops.viewpoint = 'left';
%     elc_inds = elc.global_clust_num == thisclust_global & elc.top_cons_coder &  isnan(elc.right_hemi);
%     ops.img_savename = [savedir, '/top_cons_left_clust_' num2str(iclust)];
%     ops.colormap = cons_color; 
%     plot_electrodes_on_brain_AM2(elc(elc_inds,:), ops);

    elc_inds = elc.global_clust_num == thisclust_global & elc.top_vow_coder& isnan(elc.right_hemi);
    ops.img_savename = [savedir, '/top_vow_left_clust_' num2str(iclust)];
    ops.colormap = vow_color; 
    plot_electrodes_on_brain_AM2(elc(elc_inds,:), ops);

    elc_inds = elc.global_clust_num == thisclust_global & elc.top_word_coder & isnan(elc.right_hemi);
    ops.img_savename = [savedir, '/top_syl_left_clust_' num2str(iclust)];
    ops.colormap = syl_color; 
    plot_electrodes_on_brain_AM2(elc(elc_inds,:), ops);

%         % right hemisphere
       ops.selected_hemi = 'right'; 
       ops.viewpoint = 'right';
%     elc_inds = elc.global_clust_num == thisclust_global & elc.top_cons_coder & ~isnan(elc.right_hemi);
%     ops.img_savename = [savedir, '/top_cons_right_clust_' num2str(iclust)];
%     ops.colormap = cons_color; 
%     plot_electrodes_on_brain_AM2(elc(elc_inds,:), ops);

    elc_inds = elc.global_clust_num == thisclust_global & elc.top_vow_coder & ~isnan(elc.right_hemi);
    ops.img_savename = [savedir, '/top_vow_right_clust_' num2str(iclust)];
    ops.colormap = vow_color; 
    plot_electrodes_on_brain_AM2(elc(elc_inds,:), ops);

    elc_inds = elc.global_clust_num == thisclust_global & elc.top_word_coder & ~isnan(elc.right_hemi);
    ops.img_savename = [savedir, '/top_syl_right_clust_' num2str(iclust)];
    ops.colormap = syl_color; 
    plot_electrodes_on_brain_AM2(elc(elc_inds,:), ops);

end
