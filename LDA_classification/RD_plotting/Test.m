clearvars;
load('/projectnb/busplab/Experiments/ECoG_Preprocessed_AM/electrodes.mat', 'electrodes')
load('/projectnb/busplab/Experiments/ECoG_Preprocessed_RD/filelist.mat', 'filelist');

%%% Stimuli
elecs{1} = [1; 31; 32; 35; 37; 38; 46; 47; 51; 64; 87; 103];
elecs{2} = [78; 79; 84; 123; 154; 158];
elecs{3} = [34];
elecs{4} = [58; 100; 102; 130; 199; 207];
elecs{5} = [2; 75];

cnl = do_stuff(elecs, electrodes, 2, {});
ucl = {'ESP-s'; 'ESP-r'; 'PtM-r'; 'ME-sb'; 'ME-sn'; 'AP-s'; 'No Cluster'};
for i_c = 1:length(ucl)
    cc{i_c} = sum(strcmp(ucl{i_c},cnl));
end
out = plot_on_brain(elecs, electrodes, filelist, 'stimuli_12');

%%% Subplot 1
f = figure; subplot(2,1,1);
bar(reordercats(categorical(ucl),{'ESP-s'; 'ESP-r'; 'PtM-r'; 'ME-sb'; 'ME-sn'; 'AP-s'; 'No Cluster'}),...
    cell2mat(cc))
ylabel('Stimuli')

% %%% Onset
% elecs{1} = [67; 71; 72; 103; 126; 127; 134; 136; 144; 147;];
% elecs{2} = [78; 79; 123; 154; 158];
% elecs{3} = [39];
% elecs{4} = [100; 102; 122; 130; 193; 213];
% elecs{5} = [2; 69; 75; 154; 206];
% 
% cnl = do_stuff(elecs, electrodes, 1, {});
% ucl = {'PtM-s'; 'PtM-r'; 'ME-sb'; 'ME-sn'; 'AP-r'; 'AP-s'; 'No Cluster'};
% for i_c = 1:length(ucl)
%     cc{i_c} = sum(strcmp(ucl{i_c},cnl));
% end
% plot_on_brain(elecs, electrodes, filelist, 'onset_12')
% 
% %%% Subplot 2
% figure(f); subplot(2,1,2);
% bar(reordercats(categorical(ucl),{'PtM-s'; 'PtM-r'; 'ME-sb'; 'ME-sn'; 'AP-r'; 'AP-s'; 'No Cluster'}),...
%     cell2mat(cc))
% ylabel('Onset')


function cnl = do_stuff(elecs, electrodes, elec_type, cnl)
sub_str_all = {'S357'; 'S362'; 'S369'; 'S372'; 'S376'};

for i_sub = 1:size(elecs,2)
    sub_str = sub_str_all{i_sub};
    idx = electrodes.idx_LocalProcessed(electrodes.subject == str2double(sub_str(2:end)) & ...
            electrodes.type == elec_type);
    cl = electrodes.cluster_name(electrodes.subject == str2double(sub_str(2:end)) & ...
            electrodes.type == elec_type);
    cn = cl(ismember(idx,elecs{i_sub}));
    cn(cellfun(@isempty, cn)) = {'No Cluster'};
    cnl = [cnl; cn];    
end

end

function out = plot_on_brain(elecs, electrodes, filelist, data_type)
sub_str_all = {'S357'; 'S362'; 'S369'; 'S372'; 'S376'};
base_data_path = '/projectnb/busplab/Experiments/ECoG_Preprocessed_RD/LocalEpoched/BSL1/';
addpath(genpath('/project/busplab/software/display'))

% Define colors for each subject
c = [1 0 0;... % Red
    0 1 0;... % Green
    0 0.5 1;... % Blue
    1 1 0;... % Yellow
    0 1 1]; % Teal
colors = [];
alpha = [];
shapes = [];
global_elec_idx = [];

for i_sub = 1%:size(elecs,2)
    % Calculate ei
    sub_str = sub_str_all{i_sub};
    fprintf('Processing Subject %s... ', sub_str)
    [~, preprocessed_data] = Load_Data(base_data_path, data_type, sub_str);
    data = preprocessed_data.data;
    
    s = SubjectClass(sub_str, data_type, data, filelist);
    fprintf('Finished\n');
    
    elecs_sub = elecs{i_sub};
    elecs_sub_temp = elecs_sub;
    [elecs_sub, global_elec_idx, global_elec_idx_sub] = global_elec(sub_str, electrodes, elecs_sub, global_elec_idx, data_type);
    
    ei = (abs(mean(s.W2.different.mean(elecs_sub,1000:end),2,'omitnan')) - ...
        abs(mean(s.W2.identical.mean(elecs_sub,1000:end),2,'omitnan'))) - ...
        (abs(mean(s.W1.different.mean(elecs_sub,1000:end),2,'omitnan')) - ...
        abs(mean(s.W1.identical.mean(elecs_sub,1000:end),2,'omitnan')));
      
    % Generate colors relative to how much w2 identical < w2 different
    csub = [];
    asub = [];
    ssub = [];
    for i_color = 1:length(global_elec_idx_sub)
        csub = [csub; c(i_sub,:)];
        colors = [colors; c(i_sub,:)];
        
        asub = [asub; abs(ei(i_color)/max(abs(ei)))];
        alpha = [alpha; abs(ei(i_color)/max(abs(ei)))];
        if ei(i_color) >= 0
            % Spheres
            ssub = [ssub; 's'];
            shapes = [shapes; 's'];
        else
            % Cubes
            ssub = [ssub; 'u'];
            shapes = [shapes; 'u'];
        end
    end
    fprintf('RS|RE (%s): %d | %d\n', sub_str, sum(ei >= 0), sum(ei < 0));
    % Plot only individual subject
%     plot_greenlee_electrodes_on_brain(global_elec_idx_sub, csub, asub, ssub)
    
end

out = struct(); 
out.global_elec_idx = global_elec_idx; 
out.colors = colors; 
out.alpha = alpha;
out.shapes = shapes; 

disp(shapes)
plot_greenlee_electrodes_on_brain(global_elec_idx, colors, alpha, shapes)
end

function [local_elec, global_elec_idx, global_elec_idx_sub] = global_elec(sub_str, electrodes, el, global_elec_idx, data_type)
if strcmp(data_type, 'stimuli_12')
        idx = electrodes.idx_LocalProcessed(electrodes.subject == str2double(sub_str(2:end)) & ...
            electrodes.type == 2);
        temp = idx(ismember(idx, el));
        global_elec_idx_sub = find(ismember(electrodes.idx_LocalProcessed, temp) & ...
            electrodes.subject == str2double(sub_str(2:end)) & electrodes.type==2);
    else
        idx = electrodes.idx_LocalProcessed(electrodes.subject == str2double(sub_str(2:end)) & ...
            electrodes.type == 1);
        temp = idx(ismember(idx, el));
        global_elec_idx_sub = find(ismember(electrodes.idx_LocalProcessed, temp) & ...
            electrodes.subject == str2double(sub_str(2:end)) & electrodes.type==1);
end
local_elec = el(ismember(el,temp));
global_elec_idx = [global_elec_idx;...
    global_elec_idx_sub];
end

% figure; subplot(2,1,1);
% x = reordercats(categorical({'S357', 'S362', 'S369', 'S372', 'S376'}), {'S357', 'S362', 'S369', 'S372', 'S376'});
% y = [10 -2; 4 -2; 1 0; 4 -2; 2 0];
% bar(x, y, 'stacked')
% ylabel('Stimuli')
% subplot(2,1,2)
% y = [7 -3; 3 -2; 1 0; 4 -2; 2 -3];
% bar(x, y, 'stacked')
% ylabel('Onset')