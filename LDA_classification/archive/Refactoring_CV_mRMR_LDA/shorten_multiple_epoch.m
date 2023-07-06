%{
Pass: 
N_trials dim epoch array (electrodes,3000,N_trials) with a total duration 
    of 3sec, 
an alignment condition (stimulus or onset), 
and an overall event duration (in ms) to center on the event. 

Returns a shortened epoch of duration event_duration. 
%}

function shortened_epoch = shorten_multiple_epoch(epoch_data, alignment_condition, event_duration)
num_samples = size(epoch_data, 2); 
% Calculate sample_freq in Hz. Critical that epoch_data has a duration of
% 3sec, and if not: change divisor to the epoch's duration (in seconds)
sample_freq = num_samples / 3; 

if strcmp(alignment_condition, 'stim') || strcmp(alignment_condition, 'stimulus')
    event_time = 1; % stim happens 1 second into stimulus epochs
    event_col = sample_freq * event_time; 
    num_event_samples = sample_freq * (event_duration / 1000);
    shortened_epoch = epoch_data(:,...
        event_col : event_col + num_event_samples - 1, :);    
    
elseif strcmp(alignment_condition, 'onset') 
    event_time = 2; % vox onset happens 2 seconds into onset epochs
    event_col = sample_freq * event_time; 
    num_event_samples = sample_freq * (event_duration / 1000);

    shortened_epoch = epoch_data(:,...
        (event_col - (num_event_samples / 2) + 1 :...
        event_col + (num_event_samples / 2)), :);
end


end