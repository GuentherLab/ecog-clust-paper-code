function feat_table = create_features(params, electrode_num, short_array, dynamic_timings_row)
    assert(length(size(short_array)) == 2, "epoch data for an electrode must be 2D");

    sample_freq = params.sampling_frequency;

    if strcmp(params.timing_mode, "strict")
        window_samples_ms = params.window;
        stride_samples_ms = params.stride;
    elseif strcmp(params.timing_mode, "dynamic") % && ~isnan(dynamic_timings_row)
        window_samples_ms = dynamic_timings_row.window_ms_adjusted;
        stride_samples_ms = dynamic_timings_row.stride_ms_rounded;
    end

    window_samples = sample_freq * (window_samples_ms / 1000);
    stride_samples = sample_freq * (stride_samples_ms / 1000);

    if ndims(short_array) == 2
        movmean_array = movmean(short_array, window_samples, 2, 'Endpoints', 'discard');
        stride_idxs = 1:stride_samples:size(movmean_array, 2);
        movmean_array = movmean_array(:, stride_idxs);
    end

    feat_labels_list = {};

    for window_idx = 1:size(movmean_array, 2)
        feat_labels_list{window_idx} = sprintf('e%dw%d', electrode_num, window_idx);
    end

    feat_table = array2table(movmean_array, 'VariableNames', feat_labels_list);

end
