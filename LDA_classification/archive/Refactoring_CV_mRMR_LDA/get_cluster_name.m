%{

This function returns the name of the cluster a particular electrode
belongs to, given a subject, electrode, and alignment condition.

%}

function cluster_name = get_cluster_name(sub_num_or_str, electrode_num, alignment_cond) 

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

electrode_cluster = electrode_data.cluster_name{1};
if isempty(electrode_cluster)
    cluster_name = 'NaN';
else
    cluster_name = electrode_cluster;
end

end