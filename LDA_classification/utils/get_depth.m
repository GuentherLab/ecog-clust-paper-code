%{

This function returns the depth type of a particular electrode
belongs to, given a subject, electrode, and alignment condition.

%}

function electrode_depth = get_depth(sub_num_or_str, electrode_num, alignment_cond) 

if isnumeric(sub_num_or_str)
    sub_num = sub_num_or_str;
else
    sub_str = sub_num_or_str;
    sub_num = strsplit(sub_str, 'S');
    sub_num = str2double(sub_num{2});
end

if strcmp(alignment_cond, 'stim')
    alignment_num = 2;
elseif strcmp(alignment_cond, 'onset')
    alignment_num = 1;
end

load('electrodes.mat', 'electrodes');

electrode_data = electrodes(electrodes.subject == sub_num &...
    electrodes.electrode == electrode_num &...
    electrodes.type == alignment_num, :);

electrode_depth_type = electrode_data.depth_type{1};
if isempty(electrode_depth_type)
    electrode_depth = 'NaN';
else
    electrode_depth = electrode_depth_type;
end

end