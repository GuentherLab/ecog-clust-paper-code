 %%% create tiled figures of surf plots of top coders
 % run this script after creating surf plots with surf_top_coders_per_cluster.m
 
 %% options
 
speech_feats_to_plot = {'cons', 'vow', 'syl'}; 
hems = {'left', 'right'};


% plot options
background_color = [1 1 1]; 
fig_position = [50 50 3300 730]; % [x y width height]
subplot_width = 0.15; 
subplot_height = 0.45; 
subplot_v_spacing = 0.1; % fig-relative; does not account for whitespace in images
output_resolution = 300; 

 % data processed by the script load_top_encoders.m
topelc_data_filename = [ROOT_DIR, filesep, 'projectnb/busplab/Experiments/ECoG_Preprocessed_AM/topelc_data_to_surf'];

% % % % dir to load images from
% imdir = [ROOT_DIR, filesep, 'projectnb/busplab/Experiments/ECoG_Preprocessed_AM/top_elcs_single_clust']; 
    imdir = 'C:/Users/amsme/Downloads/top_elcs_single_clust';

% % % % dir to save tiled images into
savedir = imdir; 

 

%% create plots

load(topelc_data_filename)
nfeats = length(speech_feats_to_plot);
nclusts = length(global_clust_num_list); 
nhems = length(hems);

for ifeat = [1 2 3]
    thisfeat = speech_feats_to_plot{ifeat}; 
    
    close all
    hfig = figure; 
    hfig.Color = background_color; 
    hfig.Position = fig_position; 
    
    for iclust = 1:nclusts
        for ihem = 1:nhems
            this_hem = hems{ihem};
            imfile = [imdir, '/top_', thisfeat, '_', this_hem, '_clust_', num2str(iclust), '.png'];
            im = imread(imfile);

            ypos = [subplot_height + subplot_v_spacing] * [nhems - ihem];
            xpos = [iclust - 1] / nclusts;
            %             subplot(nhems, nclusts, nclusts*[ihem-1] + iclust)
            subplot('Position',[xpos, ypos, subplot_width, subplot_height])
            image(im(200:1100, 200:1550, :))
            hax = gca;
            
            hax.XTick = [];
            hax.YTick = [];
            hax.Box = 'off';
            hax.XColor = [1 1 1];
            hax.YColor = [1 1 1];
        end 
    end
    % save the image to file
    savename = [savedir, '/tile_fig_', thisfeat];
    print(hfig, savename, '-dpng', ['-r', num2str(output_resolution)]); 
end