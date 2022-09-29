% run this section to generate/load the data 

clc
set(0,'DefaultFigureWindowStyle','docked')
addpath /project/busplab/software/ecog/scripts_LJ/Refactoring_CV_mRMR_LDA

if ~exist('p', 'var')
    p = parametersClass('/projectnb/busplab/Experiments/ECoG_fMRI_RS/Experiments/ECoG_Preprocessed_LJ/',... % path to data
        2000,...            % event duration (ms)
        50,...              % window (ms)
        25,...              % stride (ms)
        'none',...          % grouping variable
        150,...             % topN feat pooled
        30)                 % topN feat indiv
end
if ~exist('edb', 'var') 
    [~, edb] = p.generate_database; 
end
if ~exist('edb_nonan', 'var') 
    generate_edb_nonan
end

%%
% change this subset to control which electrodes get plotted:
% 2 = stim, 1 = onset
alignment = 1; 

% {'ESP-s'} = 1
% {'ESP-r'} = 2
% {'PtM-s'} = 3
% {'PtM-r'} = 4
% {'ME-sb'} = 5
% {'ME-sn'} = 6
% {'AP-r' } = 7
% {'AP-s' } = 8
global_cluster_number = 8; % clusterkey handles labeling in plot fxn

edb_nonan_subset = edb_nonan(edb_nonan.type == alignment & edb_nonan.global_clust_num == global_cluster_number,:);

% some options
ops = struct(); 
ops.force_close = 0; 
ops.clusters_as_colors = 1; 
ops.subjects_as_shapes = 1; 
ops.suppress_brain_surfs = 0; 

t = plot_electrodes_on_brain(edb_nonan_subset, ops);

% this is functioning as my shape legend for now:
unique(t)

function t = plot_electrodes_on_brain(edb_nonan_subset, ops)

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
%%% updated 2021/12/17 by Liam Jackson

addpath /project/busplab/software/ecog/util
addpath /project/busplab/software/spm12
addpath /project/busplab/software/conn
addpath /project/busplab/software/display	
addpath /project/busplab/software/display/surf

tempfile_default_name = ['mnitemp_to_delete_', datestr(datetime), '.txt']; % use new filename each time

vardefault('ops',struct);
field_default('ops','clusters_as_colors', 0)
field_default('ops','subjects_as_shapes', 0)
field_default('ops','suppress_brain_surfs', 0)
field_default('ops','base_surfaces', {'lh.pial.smoothed.surf', 'rh.pial.smoothed.surf','rh.subcortical.surf','lh.subcortical.surf'}); 
field_default('ops','tempfile', ['/projectnb2/busplab/surfshowfiles/', tempfile_default_name])
field_default('ops','img_savename',[])
field_default('ops','output_resolution',300); 
field_default('ops','force_close',0); 

clusterkey_file =  '/project/busplab/software/ecog/data/clusterkey'; 
load(clusterkey_file, 'clusterkey')
clusterkey
surfshow_resolution = 2; % 1= lowres, 2 = highres

num_electrodes = height(edb_nonan_subset);

% if toggled, map color to cluster number
if ops.clusters_as_colors
    close all % trying to make the legend with other figs open can cause problems
    color_values = edb_nonan_subset.global_clust_num; % get cluster numbers of selected electrodes
    mincolval = 1; % fix min and max vals to encompass all 8 cluster numbers
    maxcolval = 8; % fix min and max vals to encompass all 8 cluster numbers
    % make colormap for clusters_as_colors that shows each of the 8 clusters distinctly
    clusts_cmap = [ 1.0     0.5     0.0 ;...
                    0.9     0.9     0.0 ;...
                    0.1     1.0     0.0 ;...
                    0.0     0.7     1.0 ;...
                    0.0     0.0     1.0 ;...
                    0.6     0.5     1.0 ;...
                    0.8     0.0     0.8 ;...
                    1.0     0.0     0.0 ];
    field_default('ops','colormap', clusts_cmap);
    color_values = clusts_cmap(color_values, :);
else
    red_rgb = [1 0 0]; 
    color_values = repmat(red_rgb, num_electrodes, 1);
    clusts_cmap = jet(8);
    field_default('ops','colormap', clusts_cmap);
end

collapse_this = 0;
while collapse_this
% make sure there is exactly 1 color value for every electrode
% assert(num_electrodes == size(color_values, 1), 'electrode_table_indices must be same length as color_values')

% map color vals to RGB values in chosen color map
% vardefault('mincolval', min(color_values)); % set value range if not already fixed by clusters_as_colors
% vardefault('maxcolval', max(color_values)); % set value range if not already fixed by clusters_as_colors
% normvals = (color_values - mincolval) ./ (maxcolval-mincolval); % fit vals to range 0-1... baselined vals minus range
% rgbinds = 1 + round(normvals * [size(clusts_cmap,1)-1]); % match normed vals to indices in colormap
% rgbvals = ops.colormap(rgbinds,:); 
end

shapes_options = {  's',...    % - sphere (size --> side length, in mm)
                    'u',...    % - cube in MNI space (size --> diameter, in mm) || 'b',...    % - cube in image space (size --> diameter, in voxels
                    'c',...    % - crosshair in MNI space (size --> length of each axis, in mm || 'r',...    % - crosshair in image space (size --> length of each axis, in voxels
                    'i',...    % - icosahedron
                    't'};      % - tetrahedron
shapes_options_verbose = {'sphere', 'cube', 'cross', 'icosahedron', 'tetrahedron'};
                
sub_nums_vector = edb_nonan_subset.subject;
shapes = char(num_electrodes, 1);
shapes_verbose = cell(num_electrodes, 1);

if ops.subjects_as_shapes                
    unique_sub_nums = unique(sub_nums_vector); 
    for unique_sub_num_idx = 1:length(unique_sub_nums)
        shapes(sub_nums_vector == unique_sub_nums(unique_sub_num_idx)) = shapes_options{unique_sub_num_idx};
        shapes_verbose(sub_nums_vector == unique_sub_nums(unique_sub_num_idx)) = shapes_options_verbose(unique_sub_num_idx);
    end
else
    shapes(:) = shapes_options{1};
end

t = table(sub_nums_vector, shapes, shapes_verbose, edb_nonan_subset.global_clust_num,...
    'VariableNames', {'sub_num', 'shape', 'shape_long', 'global_clust_num'});

% save mni coordinates to be plotted as .txt file
mni_coords = edb_nonan_subset.mni_coord; 

% writematrix([mni_coords, zeros(size(electrode_table_indices)), shapes], ops.tempfile);
mni_coords = mat2cell(mni_coords, ones(size(mni_coords, 1),1), ones(1,size(mni_coords, 2)));
mni_coords(:,4) = mat2cell(ones(size(num_electrodes)), ones(size(num_electrodes, 1), 1), ones(1, size(num_electrodes, 2)));
mni_coords(:,5) = mat2cell(shapes, ones(size(shapes, 1), 1), ones(1, size(shapes, 2)));
writecell(mni_coords, ops.tempfile);

if ~ops.suppress_brain_surfs
    % plot electrode points with specified colors on brain surface
    surf_show('SURFACE_ADD', ops.base_surfaces,... % base surfaces
        'SURFACE_ADD', ops.tempfile, 'SURFACE_PROPERTIES','color', color_values, 'RESOLUTION', surfshow_resolution); % colored electrode locations
    delete(findobj('Type','Line')) % delete dotted outlines around patches
else
    surf_show('SURFACE_ADD', ops.tempfile, 'SURFACE_PROPERTIES','color', color_values, 'RESOLUTION', surfshow_resolution); % colored electrode locations
    delete(findobj('Type','Line')) % delete dotted outlines around patches
end

% add legend
if ops.clusters_as_colors
    add_clusters_legend(ops, clusterkey);
end

% if ops.subjects_as_shapes
%   add_shapes_legend()
% end

% if specified, save 
if ~isempty(ops.img_savename)
    surfshow_ax= gca;
    set(surfshow_ax,'Position',[0.1 0.1 .8 .8])
    tempfig= figure('Visible','off'); % Invisible figure
    temp_ax = copyobj(surfshow_ax,tempfig); % Copy the appropriate axes
    if ops.clusters_as_colors
        add_clusters_legend(ops, clusterkey)
    end
    print(tempfig, ops.img_savename, '-dpng', ['-r', num2str(ops.output_resolution)]); % save the image file
    close(tempfig)
end
    
delete(ops.tempfile) % delete temporary mni coordinates file
if ops.force_close
    close all force % close the surf_show window without requiring manual approval
end

end

% I want to combine the following into a single "generate_legend" fxn and
% accomodate for clusters_as_colors/subjects_as_shapes, haven't found a
% good way to represent the shapes in a legend so far though
%%
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

% function add_shapes_legend()
%
% end
