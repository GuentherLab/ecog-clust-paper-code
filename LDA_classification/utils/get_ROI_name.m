%{

This function returns the name of the ROI a particular electrode
belongs to, given a subject, electrode, and alignment condition.

%}

function ROI_name = get_ROI_name(sub_num_or_str, electrode_num, alignment_cond) 

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

electrode_ROI = electrode_data.roi{1};
if isempty(electrode_ROI)
    ROI_name = 'NaN';
else
    ROI_name = electrode_ROI;
end

end