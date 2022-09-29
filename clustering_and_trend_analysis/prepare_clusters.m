sub_str = {'S357', 'S362', 'S369', 'S372', 'S376'};
subj_id = 5;
subject = sub_str{subj_id};

load(strcat('/projectnb/busplab/Experiments/ECoG_Preprocessed/LocalEpoched/',subject,'/Epoch_onset_1_Hilbert_HG.mat'));
%load(strcat('/projectnb/busplab/Experiments/ECoG_Preprocessed_BG/LocalEpoched/',subject,'/Epoch_stimuli_1_LinenoiseOnly.mat'));
load(strcat('/project/busplab/software/ecog/data/',subject,'_cluster_ids.mat'));


%% Parameters
regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 10;     % maximum model order for model order estimation

samples = preprocessed_data.data;
acmaxlags = preprocessed_data.num_samples;   % maximum autocovariance lags (empty for automatic calculation)
nobs      = preprocessed_data.num_samples;   % number of observations per trial
ntrials   = size(samples, 3);

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = preprocessed_data.Fs;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

seed      = 0;      % random seed (0 for unseeded)
nvars     = 0; 

% motor_clusters = [1, 2, 5];  % stimulus clusters: ME-b, PtM-r, ME-n
motor_clusters = [1, 2, 3, 6];  % onset clusters: PtM-r, PtM-s, ME-n, ME-b


% chans = [];
% c_idx = cluster_idx{subj_id};
% for j=1:length(c_idx)
%     chans(j) = find(c_idx(j) == preprocessed_data.chan_ids);
% end
% X = samples(chans, :, :);
% for i=1:length(cluster_idx)
%     c_idx = cluster_idx{i};
%     if size(c_idx,1)==0
%         continue
%     end
%     nvars = nvars + 1;
%     chans = [];
%     for j=1:length(c_idx)
%         chans(j) = find(c_idx(j) == preprocessed_data.chan_ids);
%     end
%     X = [X; median(samples(chans, :, :), 1)];
% end
medians = median(samples, 3);
deviations = {};
c_medians = [];
c_chans = {};
key_channels = [];
key_chans = {};
chans_used = [];
for i=1:length(motor_clusters)
    c_idx = cluster_idx{motor_clusters(i)};
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
    keep = min(size(deviations{i}, 1), 4);
    chans_used = [chans_used, keep];
    key_chans{i} = chans(I(1:keep));
    key_channels = [key_channels, key_chans{i}];
end
X = samples(key_channels, :, :);

[AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
morder = moAIC;
fprintf('\nusing AIC best model order = %d\n',morder);

[A,SIG] = tsdata_to_var(X,morder,regmode);

assert(~isbad(A),'VAR estimation failed');

[G,info] = var_to_autocov(A,SIG,acmaxlags);

var_info(info,true);
F = autocov_to_pwcgc(G);
assert(~isbad(F,false),'GC calculation failed');

pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
sig  = significance(pval,alpha,mhtc);

% figure(1); clf;
% subplot(1,3,1);
% plot_pw(F);
% title('Pairwise-conditional GC');
% subplot(1,3,2);
% plot_pw(pval);
% title('p-values');
% subplot(1,3,3);
% plot_pw(sig);
% title(['Significant at p = ' num2str(alpha)])
% 
% figure(2);
% hold on
% for n=1:size(X,1)
%   plot(median(X(n,:,:),3))
% end
% legend()

results{subj_id} = struct('subject', subject, 'X', X, 'F', F, 'pval', pval, ...
'sig', sig, 'order', morder, 'info', info, 'A', A, 'SIG', SIG, 'G', G, ...
'chans_used', chans_used, 'int_key_channels', key_channels, ...
'array_channels', preprocessed_data.chan_ids(key_channels));