% PLOT_GREENLEE_ELECTRODES_ON_BRAIN: plot selected electrodes from the Greenlee/Kuzdeba dataset,
%%%%% ... with specified color values
% NOT THE SAME AS PLOT_GREENLEE_ELECTRODES_ON_BRAIN
%
% plot_greenlee_electrodes_on_brain(mni_coords_arr, color_values, ops)
%
% array of MNI coords (M x 3) for particular case
% color_values: values to be plotted (vector of same length as electrode_table_indices)...
%........ these may be any range of values; the colormap will be fit using the min and max values used here
% ops: optional struct containing fields:
%       colormap: either predefined matlab colormap or custom 256x3 RGB colormap...
%           ..... see https://www.mathworks.com/help/matlab/ref/colormap.html
%       base_surfaces: list of surfaces to plot points on top of (cell array of strings); defaults to the 4 surf_show defaults
%       tempfile: .txt file for saving mni coordinates; this file will be deleted after being loaded by surf_show
%
%%% updated 2021/8/8 by Liam Jackson

function plot_mni_coords(mni_coords_arr, cluster_color_idx, ops)

color_values = cluster_color_idx; 

vardefault('ops',struct);
field_default('ops','colormap',jet);
field_default('ops','base_surfaces', {'lh.pial.smoothed.surf', 'rh.pial.smoothed.surf','rh.subcortical.surf','lh.subcortical.surf'}); 
field_default('ops','tempfile', '/projectnb2/busplab/surfshowfiles/mnitemp.txt')

% make sure there is exactly 1 color value for every electrode
assert(size(mni_coords_arr, 1) == length(color_values), 'electrode_table_indices must be same length as color_values')

% map color vals to RGB values in chosen color map
mincolval = min(color_values);
maxcolval = max(color_values);
normvals = [color_values - mincolval] ./ maxcolval; % fit vals to range 0-1
rgbinds = 1 + round(normvals * 251); % match normed vals to indices in 1:256 colormap
rgbvals = ops.colormap(rgbinds,:); 

% save mni coordinates to be plotted as .txt file
writematrix(mni_coords_arr,ops.tempfile);

% plot elecoctrode points with specified colors on brain surface
surf_show('SURFACE_ADD', ops.base_surfaces,... % base surfaces
    'SURFACE_ADD', ops.tempfile,'SURFACE_PROPERTIES','color',rgbvals); % colored electrode locations

delete(ops.tempfile) % delete temporary mni coordinates file

