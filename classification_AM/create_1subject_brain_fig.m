 %%% create tiled figures of surf plots of top coders from subject 362
 % run this script after creating surf plots with surf_topelcs_single_subject.m
 

    close all

 %% options
 
phonunits_to_plot = {'cons', 'vow', 'syl'}; 
hems = {'left', 'right'};

% plot options
background_color = [1 1 1]; 
fig_position = [0 0 1530 1050]; % [x y width height]
% % % subplot_width = 0.15; 
% % % subplot_height = 0.45; 
% % % subplot_v_spacing = 0.1; % fig-relative; does not account for whitespace in images
% % % subplot_h_spacing = 0.1; 
output_resolution = 300; 

subject = 362; 
percent_top_elcs = 60; 
imfile_ext = '.png'; 

%%%% specify subfigure arrangement
%            phonunit, hem, viewpoint,row, column 
plottab = {'cons', 'left', 'left',   1,    1;... 
            'cons', 'right', 'right',  1,    2;...
            'cons', 'right', 'left',    1,    3;... % R hem medial
            'vowel', 'left', 'left',    2,    1;... 
            'vowel', 'right', 'right',  2,    2;...
            'vowel', 'right', 'left',    2,    3;... % R hem medial
            'syl', 'left', 'left',    3,    1;... 
            'syl', 'right', 'right',  3,    2;...
            'syl', 'right', 'left',    3,    3;... % R hem medial
            };

%% set paths, organize plot data
set_paths()

imdir = [FIG_DIR filesep 'top_elcs_single_sub'];% % % % dir to load images from
savename = [FIG_DIR filesep 'sub' num2str(subject) '_top_' num2str(percent_top_elcs) 'pct_elcs_brainplots']; %  filename to save tiled image into

plottab = cell2table(plottab, 'VariableNames', {'phonunit','hem','viewpoint','row','col'});

%% create plots
nsubplots = height(plottab); 

hfig = figure; 
hfig.Color = background_color; 
hfig.Position = fig_position; 
htl = tiledlayout(max(plottab.row),max(plottab.col));

for iplot = 1:nsubplots

    imfile = [imdir, filesep, 'sub_', num2str(subject), '_brain_top_', num2str(percent_top_elcs), 'pct_', ....
                plottab.phonunit{iplot}, '_', plottab.hem{iplot}, '_view_', plottab.viewpoint{iplot}, imfile_ext]; 
    img = imread(imfile);

    % % % % ypos = [subplot_height + subplot_v_spacing] * plottab.row(iplot) - 1;
    % % % % xpos = [subplot_width + subplot_h_spacing] * plottab.col(iplot) - 1;

    nexttile

    % image(im(200:1100, 200:1550, :))
    image(img)
    hax = gca;
    
    hax.XTick = [];
    hax.YTick = [];
    hax.Box = 'off';
    hax.XColor = [1 1 1];
    hax.YColor = [1 1 1];
    htl.Padding = 'compact';
    htl.TileSpacing = 'compact';
end

% % % % save the image to file
print(hfig, savename, '-dpng', ['-r', num2str(output_resolution)]); 