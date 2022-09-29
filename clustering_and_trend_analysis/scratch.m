subjects = {'357', '362', '369', '372', '376'};
clusters = [1, 2, 5];
% for subject=subjects
%     load(strcat('/project/busplab/software/ecog/data/S',subject{1},'_stim_cluster_ids.mat'));
%     cluster_chans = cluster_idx(clusters);
%     n_chans = 0;
%     for x=cluster_chans
%         n_chans = n_chans + size(x{1}, 2);
%     end
%     fprintf('%s: %d\n', subject{1}, n_chans);
% end
medians = median(samples, 3);
deviations = {};
c_medians = [];
c_chans = {};
key_channels = [];
key_chans = {};
for c_idx=cluster_idx(clusters) %1:length(cluster_idx)
    c_idx = c_idx{1};
    if size(c_idx,1)==0
        continue
    end
    nvars = nvars + 1;
    chans = [];
    for j=1:length(c_idx)
        chans(j) = find(c_idx(j) == preprocessed_data.chan_ids);
    end
    c_chans{i} = chans;
    c_median = median(medians(chans, :), 1);
    c_medians = [c_medians; c_median];
    deviations{i} = sum(abs(medians(chans, :) - c_median), 2);
    %[M, I] = min(deviations{i});
    [B, I] = sort(deviations{i});
    key_channels = [key_channels; chans(I)];
    key_chans{i} = chans(I);
end