addpath '/project/busplab/software/ecog/scripts_LJ/Refactoring_CV_mRMR_LDA'

p = parametersClass('/projectnb/busplab/Experiments/ECoG_fMRI_RS/Experiments/ECoG_Preprocessed_LJ/',... % path to data
    2000,...            % event duration (ms)
    50,...              % window (ms)
    25,...              % stride (ms)
    'none',...    % grouping variable
    150,...             % topN feat pooled
    30)                 % topN feat indiv

[db, edb] = p.generate_database; 

grouped_feat_struct = format_grouped_features(p, edb); 

if ~exist('all_class_sub_cv_mRMR_LDA_data', 'var')
    load(fullfile(p.grouping_path, p.grouping_variable, 'checkpoint_e2000_w50_s25.mat'), 'all_class_sub_cv_mRMR_LDA_data');
end

mdl = all_class_sub_cv_mRMR_LDA_data.word_name.S357.stim_data.word_name.stim_mdl{1}
%%
close all
clc

onset_feat_scores_all_class = table(); 
avg_scores = table();

for class_label_idx = 1:length(p.class_names)
    class_label = p.class_names{class_label_idx}; 

    feat_scores = table();

    for sub_idx = 1:length(p.sub_nums_formal)
        sub_str = p.sub_nums_formal{sub_idx};
        sub_num = p.sub_nums(sub_idx);

        topN_feat_names = horzcat(all_class_sub_cv_mRMR_LDA_data.(class_label).(sub_str).onset_data.(class_label).onset_topN_feat_names{:});
        topN_feat_scores = horzcat(all_class_sub_cv_mRMR_LDA_data.(class_label).(sub_str).onset_data.(class_label).onset_topN_scores{:});

% I'm only using the features whose predictive power (mRMR score) are 
% greater than 10% of the max value. This is slightly arbitrary and 
% needs refinement/research but it seems logical to me at least for a first pass 
        non_negligible_idxs = find(topN_feat_scores > max(max(topN_feat_scores)*0.1));
        non_negligible_idxs = 1:round((.1*length(topN_feat_scores)));
        non_negligible_idxs = 1:30;

        non_negligible_feat_scores = topN_feat_scores(non_negligible_idxs); 
        non_negligible_feat_names = topN_feat_names (non_negligible_idxs); 

        % split names into electrode id, num
        topN_electrode_names = regexprep(non_negligible_feat_names, 'w\d+', '');
        topN_electrode_nums = regexprep(topN_electrode_names, 'e', '');
        topN_electrode_nums = cell2mat(cellfun(@(x) str2double(x), topN_electrode_nums, 'UniformOutput', false));

        % split names into window id, num
        topN_window_names = regexprep(non_negligible_feat_names, 'e\d+', '');
        topN_window_nums = regexprep(topN_window_names, 'w', '');
        topN_window_nums = cell2mat(cellfun(@(x) str2double(x), topN_window_nums, 'UniformOutput', false));

        sub_feat_scores = table(categorical(topN_electrode_names', 'Ordinal', true), topN_electrode_nums', categorical(topN_window_names', 'Ordinal', true), topN_window_nums', non_negligible_feat_scores',...
            'VariableNames', {'electrode_name', 'electrode_num', 'window_name', 'window_num', 'feat_score'});

        feat_scores_temp = table(...
            {sub_feat_scores},...
            'VariableNames', {sub_str});
        feat_scores = [feat_scores, feat_scores_temp];
        feat_scores.Properties.RowNames = {class_label}; 

    end

    onset_feat_scores_all_class = [onset_feat_scores_all_class; feat_scores];
end

onset_feat_scores_all_class

for class_label_idx = 1:length(p.class_names)
    class_label = p.class_names{class_label_idx}; 

    for sub_idx = 1:length(p.sub_nums_formal)

        sub_str = p.sub_nums_formal{sub_idx};
        sub_num = p.sub_nums(sub_idx);

        sub_info = onset_feat_scores_all_class{class_label, sub_str}; 
        sub_info = sub_info{1};

        unique_electrodes = unique(sub_info.electrode_name);
        unique_windows = unique(sub_info.window_name);
        
        electrode_score_table = table();
        for electrode_idx = 1:length(unique_electrodes)
            electrode = unique_electrodes(electrode_idx); 
            electrode_avg_score = mean(sub_info.feat_score(sub_info.electrode_name == electrode));
            electrode_cumul_score = sum(sub_info.feat_score(sub_info.electrode_name == electrode));
            
            electrode_score_table_temp = table({sub_str}, electrode, electrode_avg_score, electrode_cumul_score); 
            electrode_score_table = [electrode_score_table; electrode_score_table_temp]; 
            
        end
        avg_score_table = sortrows(electrode_score_table, 'electrode_avg_score', 'ascend');
        cumul_score_table = sortrows(electrode_score_table, 'electrode_cumul_score', 'ascend');
        
        avg_electrode_sorted = avg_score_table.electrode;
        avg_score_sorted = avg_score_table.electrode_avg_score; 
        cumul_electrode_sorted = cumul_score_table.electrode; 
        cumul_score_sorted = cumul_score_table.electrode_cumul_score; 
        
        figure('IntegerHandle', 'off', 'Name', sprintf('%s, %s', class_label, sub_str));
        sgtitle(sprintf('Comparing predictive scores of top electrodes for %s, %s', sub_str, plaintext(class_label))); 
        subplot(1,3,1);
        barh(avg_score_sorted); 
        set(gca, 'YTick', 1:length(avg_electrode_sorted),...
            'YTickLabel', avg_electrode_sorted)
        title('Mean score');
        
        subplot(1,3,2);
        barh(cumul_score_sorted); 
        set(gca, 'YTick', 1:length(cumul_electrode_sorted),...
            'YTickLabel', cumul_electrode_sorted)
        title('Cumulative score');
        
        subplot(1,3,3);
        histogram(sub_info.electrode_name, 'DisplayOrder', 'ascend','Orientation','horizontal'); 
        title({'Electrode occurences in ', 'feature set'});

    end
end


%%
for sub_idx = 1:length(p.sub_nums_formal)
    sub_str = p.sub_nums_formal{sub_idx};
    sub_num = p.sub_nums(sub_idx);
    relevant_info = feat_scores.(sub_str);
    
    avg_electrode_scores = table();
    avg_window_scores = table();
    
    unique_electrodes = unique(relevant_info.electrode_name);
    unique_windows = unique(relevant_info.window_name);
    
    for electrode_idx = 1:length(unique_electrodes)
        electrode_name = unique_electrodes(electrode_idx);
        electrode_num = str2double(regexprep(electrode_name, 'e', ''));
        electrode_mean = mean(relevant_info.feat_score(strcmp(electrode_name, relevant_info.electrode_name)));
        avg_electrode_scores_temp = table(sub_num, electrode_name, electrode_num, electrode_mean, 'VariableNames', {'sub_num', 'electrode_name', 'electrode_num', 'avg_score'}); 
        avg_electrode_scores = [avg_electrode_scores; avg_electrode_scores_temp];
        
    end
    
    for window_idx = 1:length(unique_windows)
        window_name = unique_windows(window_idx);
        window_num = str2double(regexprep(window_name, 'w', ''));
        window_mean = mean(relevant_info.feat_score(strcmp(window_name, relevant_info.window_name)));
        avg_window_scores_temp = table(sub_num, window_name, window_num, window_mean, 'VariableNames', {'sub_num', 'window_name', 'window_num', 'avg_score'}); 
        avg_window_scores = [avg_window_scores; avg_window_scores_temp];

    end
end

avg_electrode_scores
avg_window_scores 




% figure();
% for sub_idx = 1:length(p.sub_nums_formal)
%     sub_str = p.sub_nums_formal{sub_idx};
%     
%     elecs = feat_scores.(sub_str).name; 
%     scores = feat_scores.(sub_str).score; 
%     
%     subplot(5, 1, sub_idx);
%     bar(scores); 
%     title(sub_str);
% end
%     


    
%%


clc 

e = {'e1', 'e256', 'e2', 'e4', 'e1'};  
w = {'w16', 'w7', 'w32', 'w50', 'w25'};  
en = [1, 256, 2, 4, 1]; 
wn = [16, 7, 32, 50, 25]; 

s = [0.01, 0.14, 0.05, 0.2, 0.1]; 

[~, ei] = sort(en); 
[~, wi] = sort(wn); 

et = table(e(ei)', en(ei)', w(ei)', wn(ei)', s(ei)', 'VariableNames', {'electrode_name', 'electrode_num', 'window_name', 'window_num', 'score'})
wt = table(e(wi)', en(wi)', w(wi)', wn(wi)', s(wi)', 'VariableNames', {'electrode_name', 'electrode_num', 'window_name', 'window_num', 'score'})


un_en = unique(en); 

%%
close all; 
clc;

sub_feat_scores
e = sub_feat_scores.electrode_name; 
en = sub_feat_scores.electrode_num; 
w = sub_feat_scores.window_name; 
wn = sub_feat_scores.window_num;

figure(1);
histogram(e, 'DisplayOrder', 'descend')

figure(2);
histogram(w, 'DisplayOrder', 'descend');

figure(3);
h2 = histogram2(en, wn); 
morebins(h2);
xlabel('en'); 
ylabel('wn');

unique_electrodes = unique(sub_feat_scores.electrode_name); 
unique_windows = unique(sub_feat_scores.window_name); 



clc
a = avg_electrode_sorted
b = cumul_electrode_sorted
c = [a, b]
a = unique(a)
b = unique(b)
c = unique(c)

%%

% for i in height edb:
%     sub num of ith row is subject to use
%     ith electrode is omitted (the ith row)
% end
% 
% rewrite all_electrodes_cv_mRMR_LDA to take in sub specific set, save accuracy only in edb compatible table
