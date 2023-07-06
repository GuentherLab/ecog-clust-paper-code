function [Word_1, Word_2, Word_Avg] = Data_Processing_Fixed_Window(data, varargin)
% Use mean of channel data as center of window and 500ms or varargin of the peak height as
% window width
if nargin == 1
    win_width = 500;
else
    win_width = varargin{1};
end

mean_word_1 = zeros(size(data,[1 2]));
mean_word_2 = zeros(size(data,[1 2]));
mean_word = zeros(size(data,[1 2]));

% win_ind_word_1 = zeros(size(data,1),2);
% win_ind_word_2 = zeros(size(data,1),2);
% win_ind = zeros(size(data,1),2);

for i = 1:size(data,1)
   % Get mean traces of each channel
   mean_word_1(i,:) = mean(data(i,:,1:2:end),3, 'omitnan');
   mean_word_2(i,:) = mean(data(i,:,2:2:end),3, 'omitnan');
   mean_word(i,:) = mean(data(i,:,:),3,'omitnan');
end

% % Get peak value and index
% [~, peak_ind_word_1] = max(mean_word_1, [], 2);
% peak_word_1 = diag(mean_word_1(:,peak_ind_word_1));
% [~, peak_ind_word_2] = max(mean_word_2, [], 2);
% peak_word_2 = diag(mean_word_2(:,peak_ind_word_2));
% [~, peak_ind] = max(mean_word, [], 2);
% peak_word = diag(mean_word(:,peak_ind));

% % Get Window indices and mean value
% win_avg_word_1 = zeros(size(data, 1), 1);
% win_avg_word_2 = zeros(size(data, 1), 1);
% win_avg = zeros(size(data, 1), 1);

% % Get average High Gamma Power from 1-3s
% hgp_win = [1000; 3000];
% hgp_word_1 = zeros(size(data, 1), 1);
% hgp_word_2 = zeros(size(data, 1), 1);
% hgp_avg = zeros(size(data, 1), 1);

% % Get Area under Curve
% win_auc_word_1 = zeros(size(data, 1), 1);
% win_auc_word_2 = zeros(size(data, 1), 1);
% win_auc = zeros(size(data, 1), 1);

% % Get Electrode Averages
% elec_avg_word_1 = zeros(size(data,1), size(data,3)/2);
% elec_avg_word_2 = zeros(size(data,1), size(data,3)/2);
% elec_avg = zeros(size(data,1), size(data,3));

% % Get Electrode Averages Over Entire Trace
% elec_avg_all_word_1 = zeros(size(data,1), size(data,3)/2);
% elec_avg_all_word_2 = zeros(size(data,1), size(data,3)/2);
% elec_avg_all = zeros(size(data,1), size(data,3));

for i = 1:size(data,1)    
    % Win Indices
%     win_ind_word_1(i,1) = max(1,peak_ind_word_1(i,:) - win_width/2);
%     win_ind_word_1(i,2) = min(size(mean_word_1,2), peak_ind_word_1(i,:) + win_width/2);
%     
%     win_ind_word_2(i,1) = max(1,peak_ind_word_2(i,:) - win_width/2);
%     win_ind_word_2(i,2) = min(size(mean_word_1,2), peak_ind_word_2(i,:) + win_width/2);
    
%     win_ind(i,1) = max(1,peak_ind(i,:) - win_width/2);
%     win_ind(i,2) = min(size(mean_word,2), peak_ind(i,:) + win_width/2);
    
%     % Win Average
%     win_avg_word_1(i) = mean(mean_word_1(i,win_ind(i,1):win_ind(i,2)), 'omitnan');
%     win_avg_word_2(i) = mean(mean_word_2(i,win_ind(i,1):win_ind(i,2)), 'omitnan');
%     win_avg(i) = mean(mean_word(i,win_ind(i,1):win_ind(i,2)), 'omitnan');
    
%     % HGP Average 1-3s
%     hgp_word_1(i) = mean(mean_word_1(i,hgp_win(1):hgp_win(2)), 'omitnan');
%     hgp_word_2(i) = mean(mean_word_2(i,hgp_win(1):hgp_win(2)), 'omitnan');
%     hgp_avg(i) = mean(mean_word(i,hgp_win(1):hgp_win(2)), 'omitnan');
    
%     % Win AUC
%     win_auc_word_1(i) = trapz(mean_word_1(i,win_ind(i,1):win_ind(i,2)));
%     win_auc_word_2(i) = trapz(mean_word_2(i,win_ind(i,1):win_ind(i,2)));
%     win_auc(i) = trapz(mean_word(i,win_ind(i,1):win_ind(i,2)));
            
%     % Elec Average
%     elec_avg_word_1(i,:) = squeeze(mean(data(i,win_ind(i,1):win_ind(i,2),1:2:end),2, 'omitnan'))';
%     elec_avg_word_2(i,:) = squeeze(mean(data(i,win_ind(i,1):win_ind(i,2),2:2:end),2, 'omitnan'))';
%     elec_avg(i,:) = squeeze(mean(data(i,win_ind(i,1):win_ind(i,2),:),2, 'omitnan'))';
    
%     % Elec Average All
%     elec_avg_all_word_1(i,:) = squeeze(mean(data(i,:,1:2:end),2, 'omitnan'))';
%     elec_avg_all_word_2(i,:) = squeeze(mean(data(i,:,2:2:end),2, 'omitnan'))';
%     elec_avg_all(i,:) = squeeze(mean(data(i,:,:),2, 'omitnan'))';
end

% % Save
Word_1.mean = mean_word_1;
% Word_1.peak = peak_word_1;
% Word_1.peak_ind = peak_ind_word_1;
% Word_1.win_ind = win_ind_word_1;
% Word_1.win_mean = win_avg_word_1;
% Word_1.hgp = hgp_word_1;
% Word_1.win_auc = win_auc_word_1;
% Word_1.elec_avg = elec_avg_word_1;
% Word_1.elec_avg_all = elec_avg_all_word_1;

Word_2.mean = mean_word_2;
% Word_2.peak = peak_word_2;
% Word_2.peak_ind = peak_ind_word_2;
% Word_2.win_ind = win_ind_word_2;
% Word_2.win_mean = win_avg_word_2;
% Word_2.hgp = hgp_word_2;
% Word_2.win_auc = win_auc_word_2;
% Word_2.elec_avg = elec_avg_word_2;
% Word_2.elec_avg_all = elec_avg_all_word_2;

Word_Avg.mean = mean_word;
% Word_Avg.peak = peak_word;
% Word_Avg.peak_ind = peak_ind;
% Word_Avg.win_ind = win_ind;
% Word_Avg.win_mean = win_avg;
% Word_Avg.hgp = hgp_avg;
% Word_Avg.win_auc = win_auc;
% Word_Avg.elec_avg = elec_avg;
% Word_Avg.elec_avg_all = elec_avg_all;
end % End Data_Processing_Numerical_Window