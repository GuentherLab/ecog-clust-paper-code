function [Word_1, Word_2, Word_Avg] = Data_Processing_Numerical_Window(data, varargin)
% Use mean of channel data as center of window and n% of the peak height as
% window width
if nargin == 1
    peak_percent = 0.7;
else
    peak_percent = varargin{1};
end

mean_word_1 = zeros(size(data,[1 2]));
mean_word_2 = zeros(size(data,[1 2]));
mean_word = zeros(size(data,[1,2]));

win_ind_word_1 = zeros(size(data,1),2);
win_ind_word_2 = zeros(size(data,1),2);
win_ind = zeros(size(data,1),2);

for i = 1:size(data,1)
   % Get mean traces of each channel
   mean_word_1(i,:) = mean(data(i,:,1:2:end),3, 'omitnan');
   mean_word_2(i,:) = mean(data(i,:,2:2:end),3, 'omitnan'); 
   mean_word(i,:) = mean(data(i,:,:),3, 'omitnan');    
end

% Get peak value and index
[~, peak_ind_word_1] = max(abs(mean_word_1), [], 2);
peak_word_1 = diag(mean_word_1(:,peak_ind_word_1));
[~, peak_ind_word_2] = max(abs(mean_word_2), [], 2);
peak_word_2 = diag(mean_word_2(:,peak_ind_word_2));
[~, peak_ind] = max(abs(mean_word), [], 2);
peak_word = diag(mean_word(:,peak_ind));

% Get Window indices and mean value
win_avg_word_1 = zeros(size(data, 1), 1);
win_avg_word_2 = zeros(size(data, 1), 1);
win_avg = zeros(size(data, 1), 1);

% Get Electrode Averages
elec_avg_word_1 = zeros(size(data,1), size(data,3)/2);
elec_avg_word_2 = zeros(size(data,1), size(data,3)/2);
elec_avg = zeros(size(data,1), size(data,3)/2);

for i = 1:size(data,1)    
    % Win Indices
%     % Word 1
%     temp = find(mean_word_1(i,1:peak_ind_word_1(i,:)) < peak_percent*peak_word_1(i,:), 1, 'last');
%     if isempty(temp)
%         win_ind_word_1(i,1) = 1;
%     else
%         win_ind_word_1(i,1) = temp;
%     end
%     
%     temp = find(mean_word_1(i,peak_ind_word_1(i,:):end) < peak_percent*peak_word_1(i,:), 1, 'first') + peak_ind_word_1(i,:) - 1;
%     if isempty(temp)
%         win_ind_word_1(i,2) = size(mean_word_1,2);
%     else
%         win_ind_word_1(i,2) = temp;
%     end
%     
%     % Word 2
%     temp = find(mean_word_2(i,1:peak_ind_word_2(i,:)) < peak_percent*peak_word_2(i,:), 1, 'last');
%     if isempty(temp)
%         win_ind_word_2(i,1) = 1;
%     else
%         win_ind_word_2(i,1) = temp;
%     end
%     
%     temp = find(mean_word_2(i,peak_ind_word_2(i,:):end) < peak_percent*peak_word_2(i,:), 1, 'first') + peak_ind_word_2(i,:) - 1;
%     if isempty(temp)
%         win_ind_word_2(i,2) = size(mean_word_2,2);
%     else
%         win_ind_word_2(i,2) = temp;
%     end
    
    % Word Avg
    temp = find(mean_word(i,1:peak_ind(i,:)) < peak_percent*peak_word(i,:), 1, 'last');
    if isempty(temp)
        win_ind(i,1) = 1;
    else
        win_ind(i,1) = temp;
    end
    
    temp = find(mean_word(i,peak_ind(i,:):end) < peak_percent*peak_word(i,:), 1, 'first') + peak_ind(i,:) - 1;
    if isempty(temp)
        win_ind(i,2) = size(mean_word,2);
    else
        win_ind(i,2) = temp;
    end
    
    % Win Average
    win_avg_word_1(i) = mean(mean_word_1(i,win_ind(i,1):win_ind(i,2)), 'omitnan');
    win_avg_word_2(i) = mean(mean_word_2(i,win_ind(i,1):win_ind(i,2)), 'omitnan');
    win_avg(i) = mean(mean_word(i,win_ind(i,1):win_ind(i,2)), 'omitnan');
    
    % Elec Average
    elec_avg_word_1(i,:) = squeeze(mean(data(i,win_ind(i,1):win_ind(i,2),1:2:end),2, 'omitnan'))';
    elec_avg_word_2(i,:) = squeeze(mean(data(i,win_ind(i,1):win_ind(i,2),2:2:end),2, 'omitnan'))';
    elec_avg(i,:) = squeeze(mean(data(i,win_ind(i,1):win_ind(i,2),2:2:end),2, 'omitnan'))';
end

% Save
Word_1.mean = mean_word_1;
Word_1.peak = peak_word_1;
Word_1.peak_ind = peak_ind_word_1;
Word_1.win_ind = win_ind_word_1;
Word_1.win_mean = win_avg_word_1;
Word_1.elec_avg = elec_avg_word_1;

Word_2.mean = mean_word_2;
Word_2.peak = peak_word_2;
Word_2.peak_ind = peak_ind_word_2;
Word_2.win_ind = win_ind_word_2;
Word_2.win_mean = win_avg_word_2;
Word_2.elec_avg = elec_avg_word_2;

Word_Avg.mean = mean_word;
Word_Avg.peak = peak_word;
Word_Avg.peak_ind = peak_ind;
Word_Avg.win_ind = win_ind;
Word_Avg.win_mean = win_avg;
Word_Avg.elec_avg = elec_avg;
end % End Data_Processing_Numerical_Window