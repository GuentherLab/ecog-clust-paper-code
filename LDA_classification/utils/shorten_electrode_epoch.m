function shortened_epoch = shorten_electrode_epoch(epoch_data, event_start_ms, ...
        event_end_ms, alignment)
    % SHORTEN_ELECTRODE_EPOCH  shortens array (epoch_data) of time-series data WRT alignment and event start/end.
    arguments
        epoch_data %(:, 3000) {mustBeNumeric}
        event_start_ms (1, 1) double {mustBeNumeric}
        event_end_ms (1, 1) double {mustBeNumeric}
        alignment {mustBeText, mustBeMember(alignment, {'stim', 'stimulus', 'onset'})}
    end

    num_samples = size(epoch_data, 2);
    assert(num_samples == 3000, "Epoch data must be full length (3sec of samples @ 1kHz = 3000samples)");
    % Calculate sample_freq in Hz. Critical that epoch_data has a duration of
    % 3sec, and if not: change divisor to the epoch's duration (in seconds)
    sample_freq = parametersClass.sampling_frequency;

    if strcmp(alignment, 'stim') || strcmp(alignment, 'stimulus')
        event_time = 1; % stim happens 1 second into stimulus epochs
        event_col = sample_freq * event_time;
    elseif strcmp(alignment, 'onset')
        event_time = 2; % vox onset happens 2 seconds into onset epochs
        event_col = sample_freq * event_time;
    end

    event_start_sample = sample_freq * (event_start_ms / 1000); % Hz * (ms/1000) = Hz * sec
    event_start = event_col + event_start_sample +1;

    event_end_sample = sample_freq * (event_end_ms / 1000);
    event_end = event_col + event_end_ms - 1;

    shortened_epoch = epoch_data(:, ...
        event_start:event_end);

end
