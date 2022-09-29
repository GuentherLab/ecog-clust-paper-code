sub_str = {'357', '362', '369', '372', '376'};
subj_id = 1;
subject = 'S357';

load(strcat('/projectnb/busplab/Experiments/ECoG_Preprocessed_BG/LocalEpoched/',subject,'/Epoch_stimuli_1_LinenoiseOnly.mat'));
load(strcat('/project/busplab/software/ecog/data/',subject,'_stim_cluster_ids.mat'));


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
fres      = 1000;     % frequency resolution (empty for automatic calculation)

seed      = 0;      % random seed (0 for unseeded)
nvars     = 0; 

X = [];
for i=1:length(cluster_idx)
    c_idx = cluster_idx{i};
    if size(c_idx,1)==0
        continue
    end
    nvars = nvars + 1;
    chans = [];
    for j=1:length(c_idx)
        chans(j) = find(c_idx(j) == preprocessed_data.chan_ids);
    end
    X = [X; median(samples(chans, :, :), 1)];
end


n = size(X,1);
step = 10;
window = 100;
n_samples = size(X,2);
n_frames = floor((n_samples - window) / step) + 1;
St = zeros(n, n, fres+1, n_frames);
Zt = zeros(n, n, fres+1, n_frames);
for k = 1:step:n_samples-window
    [Z,S] = tsdata_to_erc(X(:, k:window+k, :), momax, fres);
    Zt(:, :, :, k) = Z;
    St(:, :, :, k) = S;
end