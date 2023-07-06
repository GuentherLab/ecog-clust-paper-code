% PLOT_GREENLEE_ELECTRODES_ON_BRAIN: plot selected electrodes from the Greenlee/Kuzdeba dataset,
%%%%% ... with specified color values
%
% plot_greenlee_electrodes_on_brain(electrode_table_indices, color_values, ops)
%
% electrode_table_indices: list of indices to plot from /projectnb2/busplab/Experiments/ECoG_Preprocessed_AM/electrodes.mat
% color_values: values to be plotted (vector of same length as electrode_table_indices)...
%........ these may be any range of values; the colormap will be fit using the min and max values used here
% ops: optional struct containing fields:
%       colormap: either predefined matlab colormap or custom 256x3 RGB colormap...
%           ..... see https://www.mathworks.com/help/matlab/ref/colormap.html
%       clusters_as_colors: if true, assign colors based on cluster, using electrodes.global_clust_num... 'color_values' will be ignored
%          ... cluster colors will be consistent within a colormap and across alignment conditions,
%          ... based on Kuzdeba thesis table 1 start time
%       base_surfaces: list of surfaces to plot points on top of (cell array of strings); defaults to the 4 surf_show defaults
%       tempfile: .txt file for saving mni coordinates; this file will be deleted after being loaded by surf_show
%       img_savename: if not empty, the brain image will be saved as a .png with this filename
%       output_resolution: if img_savename is not empty, the saved image will be saved as this resolution in dots per inch (default 300)
%       force_close: if true, immediately close the surf_show window with 'close all force' (useful running this func in a batch to save figures)
%
% LJ Update: seeks to create functionality to just send in a list of
% electrodes to plot

%%% updated 2021/12/15 by Liam Jackson

function plot_greenlee_electrodes_on_brain(edb_nonan, electrode_table_indices, color_values, ops)

addpath /project/busplab/software/ecog/util/
addpath /project/busplab/software/spm12
addpath /project/busplab/software/conn
addpath /project/busplab/software/display	
addpath /project/busplab/software/display/surf

tempfile_default_name = ['mnitemp_to_delete_', datestr(datetime), '.txt']; % use new filename each time

electrodes_table_file = '/projectnb2/busplab/Experiments/ECoG_Preprocessed_AM/electrodes.mat'; 
clusterkey_file =  '/project/busplab/software/ecog/data/clusterkey'; 
vardefault('ops',struct);
field_default('ops','clusters_as_colors', 0)
field_default('ops','base_surfaces', {'lh.pial.smoothed.surf', 'rh.pial.smoothed.surf','rh.subcortical.surf','lh.subcortical.surf'}); 
field_default('ops','tempfile', ['/projectnb2/busplab/surfshowfiles/', tempfile_default_name])
field_default('ops','img_savename',[])
field_default('ops','output_resolution',300); 
field_default('ops','force_close',0);
field_default('ops','img_title','');

load(electrodes_table_file, 'electrodes') % electrodes data table
electrodes = edb_nonan;
load(clusterkey_file, 'clusterkey')

surfshow_resolution = 2; % 1= lowres, 2 = highres

% if toggled, map color to cluster number
if ops.clusters_as_colors
    close all % trying to make the legend with other figs open can cause problems
    non_nan_inds = find(~isnan(electrodes.global_clust_num(electrode_table_indices)))
    electrode_table_indices = electrode_table_indices(non_nan_inds); % get rid of electrodes with no cluster label
    color_values = electrodes.global_clust_num(electrode_table_indices); % get cluster numbers of selected electrodes
    mincolval = 1; % fix min and max vals to encompass all 8 cluster numbers
    maxcolval = 8; % fix min and max vals to encompass all 8 cluster numbers
    % make colormap for clusters_as_colors that shows each of the 8 clusters distinctly
    clusts_cmap = [   1, 0.5,   0;...
                    0.9, 0.9,   0;...
                    0.1,   1,   0;...
                      0, 0.7,   1;...
                      0,   0,   1;...
                    0.6, 0.5,   1;...
                    0.8,   0, 0.8;...
                      1,   0,   0];
                  
    field_default('ops','colormap', clusts_cmap);
else
    field_default('ops','colormap', jet)
end

% make sure there is exactly 1 color value for every electrode
assert(length(electrode_table_indices) == length(color_values), 'electrode_table_indices must be same length as color_values')

% map color vals to RGB values in chosen color map
vardefault('mincolval', min(color_values)); % set value range if not already fixed by clusters_as_colors
vardefault('maxcolval', max(color_values)); % set value range if not already fixed by clusters_as_colors
normvals = (color_values - mincolval) ./ (maxcolval-mincolval); % fit vals to range 0-1... baselined vals minus range
rgbinds = 1 + round(normvals * (size(clusts_cmap,1)-1)); % match normed vals to indices in colormap
rgbvals = ops.colormap(rgbinds,:); 

% save mni coordinates to be plotted as .txt file
mni_coords = electrodes.mni_coord(electrode_table_indices,:); 
writematrix(mni_coords,ops.tempfile);

% plot elecoctrode points with specified colors on brain surface
% surf_show('SURFACE_ADD', ops.base_surfaces, 'SURFACE_PROPERTIES', 'transparency', .5, 'RESOLUTION', surfshow_resolution); % base surfaces
subsurf_transparency = 0.15;
surf_show('SURFACE_ADD', 'lh.pial.smoothed.surf', 'SURFACE_PROPERTIES', 'transparency', subsurf_transparency, 'RESOLUTION', surfshow_resolution,...
    'SURFACE_ADD', 'lh.subcortical.surf', 'SURFACE_PROPERTIES', 'transparency', subsurf_transparency, 'RESOLUTION', surfshow_resolution,...
    'SURFACE_ADD', ops.tempfile, {}, true, true, 'SURFACE_PROPERTIES','color',rgbvals, 'RESOLUTION', surfshow_resolution); % colored electrode locations
delete(findobj('Type','Line')) % delete dotted outlines around patches

% add legend
if ops.clusters_as_colors
    add_clusters_legend(ops, clusterkey)
end

% if specified, save 
if ~isempty(ops.img_savename)
    surfshow_ax= gca;
    set(surfshow_ax,'Position',[0.1 0.1 .8 .8])
    title(ops.img_title);
    tempfig= figure('Visible','off'); % Invisible figure
    temp_ax = copyobj(surfshow_ax,tempfig); % Copy the appropriate axes
    if ops.clusters_as_colors
        add_clusters_legend(ops, clusterkey)
    end
    print(tempfig, ops.img_savename, '-dpng', ['-r', num2str(ops.output_resolution)]); % save the image file
    close(tempfig)
end
    
if ops.force_close
    close all force % close the surf_show window without requiring manual approval
end

delete(ops.tempfile) % delete temporary mni coordinates file
end

function add_clusters_legend(ops, clusterkey)
    allpatches = table(findobj('Type','Patch'), 'VariableNames',{'patch'});
    npatches = height(allpatches);
    allpatches.global_clust_num = NaN(npatches,1); 
    % in surf show, skip first 2 and last 2 Patch objects, which are something other than the points we added
%     for ipatch = 3:npatches-2
    for ipatch = 1:npatches
        thiscol = allpatches.patch(ipatch).FaceVertexCData(1,:); % color of this patch 
        thisclust = find(all(thiscol == ops.colormap,2)); 
        if ~isempty(thisclust)
            allpatches.global_clust_num(ipatch) = thisclust; % find which clustnum produced this patch color
        end
    end
    [clustnums, inds4leg] = unique(allpatches.global_clust_num);
    inds4leg = inds4leg(~isnan(clustnums)); % exclude patches which are not points we plotted
    clustnums = clustnums(~isnan(clustnums)); 
    nclusts = length(clustnums);
    leglabels = cell(size(clustnums)); 
    for iclust = 1:nclusts % look up the cluster names from clusterkey
        clusterkey_matchrow = find(clusterkey.global_clust_num==clustnums(iclust),1);
        leglabels{iclust} = clusterkey.name{clusterkey_matchrow}; 
    end
    legend(allpatches.patch(inds4leg),leglabels)
end

