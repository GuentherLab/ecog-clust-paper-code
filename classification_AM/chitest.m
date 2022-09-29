function [chisq_pval] = chitest(top_elcs_per_cluster, total_elcs_per_cluster, top_elcs_proportion, expected_elcs)

% updated by Andrew Meier 2022/7/30
%%% leave top_elcs_proportion as [] if predicted_elcs vector is specified

assert( isempty(top_elcs_proportion) | ~exist('expected_elcs','var'),...
    'leave top_elcs_proportion as [] if expected_elcs vector is specified')

% default for expected_elcs: top electrodes if they were proportionately distributed across clusters
if ~exist('expected_elcs','var')
    expected_elcs = top_elcs_proportion * total_elcs_per_cluster; 
end
    
noutcomes = length(top_elcs_per_cluster);
degfr = noutcomes-1;
devstat = (top_elcs_per_cluster-expected_elcs).^2 ./ expected_elcs; % calculate deviations from expected for each outcome
chistat = sum(devstat); % chi2 test statistic
chisq_pval = chi2cdf(chistat,degfr,'upper'); % get pval





%%%%%% note: this function appears mostly redundant with a matlab built-in:
% [h p stats] = chi2gof([1:length(top_elcs_per_cluster)]', 'Frequency',top_elcs_per_cluster, 'Expected',expected_elcs, 'Emin',0)