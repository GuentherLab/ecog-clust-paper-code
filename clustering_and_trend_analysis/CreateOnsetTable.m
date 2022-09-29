clear all; clc

sub_str = 'S376'; % {'S352', 'S357', 'S362', 'S369', 'S372', 'S376'} 
sess_str = ''; % {'', '2'}

% trig_Fs = 24400; % trigger sampled at 24kHz
switch sub_str
    case 'S352'
        trig_Fs = 24000;
    case 'S357'
        trig_Fs = 24000;
    case 'S362'
        trig_Fs = 24000;
    case 'S369'
        trig_Fs = 16000;
    case 'S372'
        trig_Fs = 16000;
    case 'S376'
        trig_Fs = 16000;
end

base_path = error('Enter data directory path here');

try
    load(fullfile(base_path, 'ProcessedDataset', [sub_str add_str 'OnsetTime.mat']));
catch
    load(fullfile(base_path, 'LocalProcessed', sub_str, [sub_str add_str 'OnsetTime.mat']));
end

OnsetTable = cell(0);
for i=1:length(stimuliOnset)

    stimuli_t = stimuliOnset(i).VisualOnset;
    onset_t = stimuliOnset(i).SpeechOnset;
    
    % Convert to ms (speech and trigger sampled at 24.4kHz)
    stimuli_t = round(1000 * stimuli_t/trig_Fs);
    onset_t = round(1000 * onset_t/trig_Fs);
    
    % Line up speech with triggers/stimuli
    matched_onset_t = nan(size(stimuli_t));
    for k=1:length(stimuli_t)-1
        match_ind = find(onset_t > stimuli_t(k) & onset_t < stimuli_t(k+1));
        if isempty(match_ind) || length(match_ind) > 1
            matched_onset_t(k) = nan;
        else
            matched_onset_t(k) = onset_t(match_ind);
        end
    end
    
    onset_table = [stimuli_t(1:2:end)' stimuli_t(2:2:end)' ...
        matched_onset_t(1:2:end)' matched_onset_t(2:2:end)'];
    
    bad_trial_1 = find(isnan(onset_table(:,3)));
    bad_trial_2 = find(isnan(onset_table(:,4)));
    
    % First speech is bad, entire trial bad
    onset_table(bad_trial_1, 1) = nan;
    onset_table(bad_trial_1, 2) = nan;
    % Second speech is bad, keep first part (no repetition though)
    onset_table(bad_trial_2, 2) = nan;
    
    OnsetTable{i} = onset_table;
end


save_file = fullfile(base_path, 'LocalProcessed', sub_str, ['LocalOnsetTable' sess_str '.mat']);
save(save_file, 'OnsetTable');
    
