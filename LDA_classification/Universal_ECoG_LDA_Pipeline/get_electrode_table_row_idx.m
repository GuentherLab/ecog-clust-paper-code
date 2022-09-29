%{

This function returns the row number from the electrode electrode table
that an electrode belongs to, given a subject, electrode id, and alignment 
condition.

%}

function elect_table_row_idx = get_electrode_table_row_idx(sub_num_or_str, electrode_num, alignment_cond) 

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

elect_table_row_idx = find(electrodes.subject == sub_num &...
    electrodes.electrode == electrode_num &...
    electrodes.type == alignment_num);

end