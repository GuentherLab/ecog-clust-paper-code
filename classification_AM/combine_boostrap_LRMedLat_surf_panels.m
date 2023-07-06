 %%%% combine the surf brainplots created by areal_phonunit_preference_1elc_included.m into a single tiled figure
 % this plot will include L lat, L med, R lat, R med views
%%%% before running this script, run areal_phonunit_preference_1elc_included.m and move outputs to the directory specified below

 clear
 set_paths()

 %% options
 %%%%%% output figure will be saved into bootstrap_results_dir



 %%%%% results from 2-second window analysis, with Bonferroni correction, with  ops.phon_pref_ind_mode = 'greater_than_one'
 ops.bootstrap_results_dir = [DIR_LDA_2SEC_WINDOW filesep 'PPI greater than 1, bonf cor'];

  %%%%%% results from 2-second window analysis, without Bonferroni correction, with  ops.phon_pref_ind_mode = 'greater_than_one'
% ops.bootstrap_results_dir = [DIR_LDA_2SEC_WINDOW filesep 'PPI greater than 1, no bonf cor'];



% plot options
ops.background_color = [1 1 1]; 
ops.fig_position = [50 50 1370 950]; % [x y width height]
ops.output_resolution = 300; 
ops.edge_padding = 'Tight'; % space between figure borders and panels

ops.tile_spacing = 'Tight'; % space between panels.... 'Loose','Compact','Tight','None'
% ops.tile_spacing = 'Compact'; % space between panels.... 'Loose','Compact','Tight','None'

hem_file_ending_strs = {'left_hem_left_view.png';...
                        'right_hem_right_view.png';...
                        'left_hem_right_view.png';...
                        'right_hem_left_view.png'};


 %% load images, make figure
 % data_filename = 'single_elc_included_trialwise_accuracy.mat';
 % load([bootstrap_results_dir, filesep, data_filename]) % maybe don't need to load this data

close all

file_leading_string = [ops.bootstrap_results_dir filesep 'phon_pref_brainplot_all_subs_']; 

hfig = figure; 
hfig.Color = ops.background_color; 
hfig.Position = ops.fig_position; 

htile = tiledlayout(2,2);
htile.Padding = ops.edge_padding;
htile.TileSpacing = ops.tile_spacing; 

for ifile = 1:size(hem_file_ending_strs,1)
    [hax, hplotimg] = plot_and_format_img([file_leading_string, hem_file_ending_strs{ifile}], ops);
end






%% save figure
savename = [ops.bootstrap_results_dir filesep 'phon_pref_brainplot_all_subs_all_views']; 
print(hfig, savename, '-dpng', ['-r', num2str(ops.output_resolution)]); 

%%
function [hax, hplotimg] = plot_and_format_img(hem_file_str, ops)
    nexttile
    img = imread(hem_file_str);
    hplotimg = image(img);

    hax = gca;
    hax.XTick = [];
    hax.YTick = [];
    hax.Box = 'off';
    hax.XColor = [1 1 1];
    hax.YColor = [1 1 1];
end
