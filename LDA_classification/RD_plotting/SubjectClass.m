classdef SubjectClass
    
    %% Public Properties
    properties
       ID
       Type
       Data
       
       Trials
       Labels
       Conditions
       
       Depth_Type
       Cluster
       
       W1
       W2
       WAvg
       
       EI
       EI_SD
       P
       P_Test
       W_Test
       AE       
    end % End Properties
    
    %% Public Methods    
    methods
        
        % Constructor
        function obj = SubjectClass(sub_str, data_type, data, filelist)
            obj.ID = sub_str;
            obj.Type = data_type;
            obj.Data = data;  
            
            [obj.Trials, obj.Labels, obj.Conditions] = Find_Trial_Details(obj, filelist);
            [obj.Depth_Type, obj.Cluster] = Find_Electrode_Details(obj);
            
            [obj.W1, obj.W2, obj.WAvg] = Generate_W(obj);
            
            [obj.EI, obj.EI_SD] = Calculate_EI(obj);
%             obj.P = Permutation_Test(obj); %Permutation_Test(obj); % Calculate_Electrode_P(obj);
%             obj.P_Test = Permutation_Test(obj);
%             obj.W_Test = Wilcoxon_Test(obj);
            obj.AE = Active_Electrodes(obj);
        end
        
    end % End Public Methods
    
    %% Private Methods
    methods (Access = private)
        
        %% Get details on subject trial data 
        function [trials, labels, conditions] = Find_Trial_Details(obj, filelist)
            filelist_sub = filelist(find(filelist.subject==str2double(obj.ID(2:end))),:);
            trials = filelist_sub.trials{1,1};
            labels = trials.pair_comparison;                
            conditions = sort(unique(labels));
        end
        
        %% Get details on subject electrode data
        function [depth_type, cluster] = Find_Electrode_Details(obj)
            load('/projectnb/busplab/Experiments/ECoG_Preprocessed_AM/electrodes.mat', 'electrodes')
            idx = electrodes.idx_LocalProcessed(electrodes.subject == str2double(obj.ID(2:end)) & electrodes.type == strcmp(obj.Type, 'stimuli_12')+1);
                        
            temp = electrodes.cluster_name(electrodes.subject == str2double(obj.ID(2:end)) & electrodes.type == strcmp(obj.Type, 'stimuli_12')+1);
            cluster = cell(size(obj.Data,1), 1);
            cluster(idx) = temp;
            cluster(cellfun(@isempty, cluster)) = {'No_Cluster'};
                        
            temp = electrodes.depth_type(electrodes.subject == str2double(obj.ID(2:end)) & electrodes.type == strcmp(obj.Type, 'stimuli_12')+1);
            depth_type = cell(size(obj.Data,1), 1);
            depth_type(idx) = temp;
            depth_type(cellfun(@isempty, depth_type)) = {'No_Depth_Type'};
        end
        
        %% Generate W1 and W2
        function [w1, w2, wavg] = Generate_W(obj)
            [w1.all, w2.all, wavg.all] = Data_Processing_Fixed_Window(obj.Data);                
            
            for i = 1:size(obj.Conditions,1)
                [w1.(obj.Conditions{i}), w2.(obj.Conditions{i}), wavg.(obj.Conditions{i})] ...
                    = Data_Processing_Fixed_Window(obj.Data(:,:,find(strcmp(obj.Trials.pair_comparison,obj.Conditions{i}))));
            end            
        end
        
        %% Calculate EI
        function [ei, ei_sd] = Calculate_EI(obj)
            temp = mean(obj.Data(:,1001:3000,:), 2, 'omitnan');
            ei.all = mean(temp(:,:,2:2:end) - temp(:,:,1:2:end), 3, 'omitnan');
            ei_sd.all = std(temp(:,:,2:2:end) - temp(:,:,1:2:end), 0, 3, 'omitnan');
            
            for i_cond = 1:size(obj.Conditions,1)
                temp = mean(obj.Data(:,1001:end,strcmp(obj.Labels, obj.Conditions{i_cond})), 2, 'omitnan');
                ei.(obj.Conditions{i_cond}) = mean(temp(:,:,2:2:end) - temp(:,:,1:2:end), 3, 'omitnan');
                ei_sd.(obj.Conditions{i_cond}) = std(temp(:,:,2:2:end) - temp(:,:,1:2:end), 0, 3, 'omitnan');
            end
        end
        
        %% Identify Significant Electrodes
        function p = Calculate_Electrode_P(obj)
           elec = cell(size(obj.Conditions));
           p = zeros(size(obj.Data, 1), 1);
           
           for i_cond = 1:size(obj.Conditions, 1)
              temp = mean(obj.Data(:,1001:end,strcmp(obj.Labels, obj.Conditions{i_cond})), 2, 'omitnan');
              elec{i_cond} = squeeze(temp); %squeeze(temp(:,:,2:2:end) - temp(:,:,1:2:end));
           end
           
           for i_elec = 1:size(obj.Data ,1)
               conds = [1; 3]; % Only look at Different and Identical
               temp = zeros(size(elec{1}, 2)/2, size(conds, 1));
               
               % On all conditions (1:size(...)) or
               % on just Identical and different (1:2:size(...))
               % Word 1 not significant?
               for i_cond = 1:size(conds, 1)
                   % If there are a different number of trial pairs
                   try
                      temp(:, i_cond) = elec{conds(i_cond)}(i_elec,1:2:end)';
                   catch
                      min_elecs = min(size(temp,1), size(elec{conds(i_cond)}(i_elec,1:2:end)', 1));
                      temp(1:min_elecs, i_cond) = elec{conds(i_cond)}(i_elec,1:min_elecs)';
                   end
               end
               [p_temp, ~, ~] = kruskalwallis(temp, [], 'off');
               
               % Word 2 significant?
               if p_temp > 0.05 % If no significant difference in w1
                  for i_cond = 1:size(conds, 1)
                      % If there are a different number of trial pairs
                      try
                         temp(:, i_cond) = elec{conds(i_cond)}(i_elec,2:2:end)';
                      catch
                         min_elecs = min(size(temp,1), size(elec{conds(i_cond)}(i_elec,2:2:end)', 1));
                         temp(1:min_elecs, i_cond) = elec{conds(i_cond)}(i_elec,1:min_elecs)';
                      end
                  end
                   
                  [p(i_elec), ~, ~] = kruskalwallis(temp, [], 'off'); 
               end           
           end
        end
        
        %% Permutation Test
        function p = Permutation_Test(obj)
            p = zeros(size(obj.Data,1),1);
            
            % Variable parameters
            conds = [1; 3]; % Conditions being looked at (1-Different, 2-Flipped, 3-Identical)
            nperm = 1e4;
            
            % For each electrode
            for i_elec = 1:size(obj.Data,1)
                % Separate data into conditions (153x3000x216 -> 72x2000 for each cell index)
                temp = cell(size(conds,1), 1);
                for i_cond = 1:size(conds,1)
                    temp{i_cond} = squeeze(obj.Data(i_elec,1001:end,strcmp(obj.Labels, obj.Conditions{conds(i_cond)})))';
                end
                
                % Ensure conditions have the same number of trials
                M = min(size(temp{1},1), size(temp{2},1));
                if size(temp{1},1) ~= size(temp{2},1)                   
                   temp{1} = temp{1}(1:M,:);
                   temp{2} = temp{2}(1:M,:);
                end
                
                % Word 1
                elec_data = [temp{1}(1:2:end,:); temp{2}(1:2:end,:)];

                % Permutation Test
                fun = @(elec_data) norm(mean(elec_data(1:M/2,:),1,'omitnan') - mean(elec_data(M/2+1:end,:),1,'omitnan'));
                d1 = fun(elec_data);
                d0 = arrayfun(@(n) fun(elec_data(randperm(size(elec_data,1)),:)), 1:nperm);
                p_temp = mean(d0>d1,'omitnan');
                
                % Word 2
                if p_temp > 0.05
                    elec_data = [temp{1}(2:2:end,:); temp{2}(2:2:end,:)];

                    % Permutation Test
                    fun = @(elec_data) norm(mean(elec_data(1:M/2,:),1,'omitnan') - mean(elec_data(M/2+1:end,:),1,'omitnan'));
                    d1 = fun(elec_data);
                    d0 = arrayfun(@(n) fun(elec_data(randperm(size(elec_data,1)),:)), 1:nperm);
                    p(i_elec) = mean(d0>d1,'omitnan');
                end
            end                        
        end
        
        %% Wilcoxon Signed Rank Test
        function p = Wilcoxon_Test(obj)            
            p = zeros(size(obj.Data,1),1);
            
            % Variable parameters
            conds = [1; 3]; % Conditions being looked at (1-Different, 2-Flipped, 3-Identical)
            tstart = 1001;
            
            % For each electrode
            for i_elec = 1:size(obj.Data,1)
                % Separate data into conditions (153x3000x216 -> 72x2000 for each cell index)
                temp = cell(size(conds,1), 1);
                for i_cond = 1:size(conds,1)
                    temp{i_cond} = squeeze(obj.Data(i_elec,tstart:end,strcmp(obj.Labels, obj.Conditions{conds(i_cond)})))';
                end
                
                % Ensure conditions have the same number of trials
                M = min(size(temp{1},1), size(temp{2},1));
                if size(temp{1},1) ~= size(temp{2},1)                   
                   temp{1} = temp{1}(1:M,:);
                   temp{2} = temp{2}(1:M,:);
                end
                
                % Take Wilcoxon Signed Rank test for difference between
                % conditions (null hypothesis of zero mean diff)
                elec_diff = temp{1} - temp{2};                
                p(i_elec) = signrank(mean(elec_diff,1,'omitnan'));
%                 p_temp = zeros(size(elec_diff,1),1);
%                 for i = 1:size(elec_diff,1)
%                     if sum(isnan(elec_diff(i,:))) == size(elec_diff,2)
%                         p_temp(i) = NaN;
%                     else
%                         p_temp(i) = signrank(elec_diff(i,:));
%                     end
%                 end
%                 
%                 % Average p_values for each trial
%                 p(i_elec) = mean(p_temp, 'omitnan');
                
                if p(i_elec) < 0.05
%                    figure; hold on;
%                    histogram(p_temp, max(floor(length(p_temp)/5), 6))
%                    xline(0, '--r');
%                    yline(0, '-k');
%                    title([obj.ID ': ' num2str(i_elec) newline...
%                        ' Distribution of p-values Across Trials'])
                end
            end
            
            
            
            
        end
        
        %% Identify Electrodes showing any Task Related Activity
        function ae = Active_Electrodes(obj)
            win = 250;
            start = 2000;
            stop = 3000;
            threshold = 1;
            
            nDur = 100;

            % Truncate to window
            % Stimuli Case
            if strcmp(obj.Type, 'stimuli_12')
                data = obj.Data(:,start:stop,:);
            
            % Onset case
            else
                if 0 % Use 2000-3000ms window in Onset
                    data = obj.Data(:,start:stop,:);
                else % Use same time window as Stimuli           
                    load(['/projectnb2/busplab/Experiments/ECoG_Preprocessed_RD/LocalProcessed/' obj.ID '/LocalOnsetTable.mat'])
                    onset_table = OnsetTable{1,1};
                    
                    vo_line = 2000;
                    so_line_1 = vo_line - (onset_table(:,3)-onset_table(:,1));
                    so_line_2 = vo_line - (onset_table(:,4)-onset_table(:,2));
                    
                    for i = 1:length(so_line_1)
                       start(2*i - 1) = so_line_1(i);
                       start(2*i) = so_line_2(i);                        
                    end
                    
                    stop = start+1000;
                    try
                        data = obj.Data(:,start:min(stop,3000),:); 
                    catch
                        pause
                    end
                end
            end % End truncate window

            % Split dataset by conditions and word, then average
            for i_cond = 1:size(obj.Conditions,1)% different, flipped, identical
                idx = strcmp(obj.Trials.pair_comparison, obj.Conditions{i_cond});
                temp = data(:, :, idx);

                w1(:, :, i_cond) = mean(temp(:, :, 1:2:end), 3, 'omitnan');
                w2(:, :, i_cond) = mean(temp(:, :, 2:2:end), 3, 'omitnan'); 
            end

            % Process
            chans = cell(size(obj.Conditions,1),1);
            temp = [];
            for i_cond = 1:size(obj.Conditions,1)
                % word 1
                for i_win = 1:size(w1,2)-win
                    temp(:,i_win) = abs(mean(w1(:,i_win+win,i_cond),2,'omitnan')) > threshold;
                end
                chans1 = sum(temp,2) > nDur;

                % word 2
                for i_win = 1:size(w2,2)-win
                    temp(:,i_win) = abs(mean(w2(:,i_win+win,i_cond),2,'omitnan')) > threshold;
                end
                data_chans2 = find(sum(temp,2) > nDur);
                chans2 = sum(temp,2) > nDur;

                % channels in either w1 and w2
                chans{i_cond} = chans1 | chans2;
            end

            % channels in both w1 and w2 across any condition
            ae = chans{1} | chans{2} | chans{3};  
        end
    end % End Private Methods
end % End Classdef