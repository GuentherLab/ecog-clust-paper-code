% window state: surf_show line 90
close all force; 

class_label = p.class_names{3};


left_stim_files = dir(strcat('*stim*', class_label, '*.png'))

left_idx = [1 4 2 3];

figure();
tiles = tiledlayout(2,4);
tiles.TileSpacing = 'tight'; 
tiles.Padding = 'tight';

for file = left_idx
nexttile
imshow(left_stim_files(file).name)
end

onset_files = dir(strcat('*onset*', class_label, '*.png'))

left_idx = [1 4 2 3];

% figure();
% tiles = tiledlayout(2,2);
% tiles.TileSpacing = 'tight'; 
% tiles.Padding = 'tight';

for file = left_idx
nexttile
imshow(onset_files(file).name)
end


%%
close all force; 

files = {   'top_20_pct_stim_consonants_name_left_hemi_left_view.png',
            'top_20_pct_stim_consonants_name_right_hemi_right_view.png',
            'top_20_pct_stim_vowel_name_left_hemi_left_view.png',
            'top_20_pct_stim_vowel_name_right_hemi_right_view.png',
            'top_20_pct_stim_word_name_left_hemi_left_view.png',
            'top_20_pct_stim_word_name_right_hemi_right_view.png',
            'top_20_pct_onset_consonants_name_left_hemi_left_view.png',
            'top_20_pct_onset_consonants_name_right_hemi_right_view.png',
            'top_20_pct_onset_vowel_name_left_hemi_left_view.png',
            'top_20_pct_onset_vowel_name_right_hemi_right_view.png',
            'top_20_pct_onset_word_name_left_hemi_left_view.png',
            'top_20_pct_onset_word_name_right_hemi_right_view.png'};
figure();
tiles = tiledlayout(2,6);
tiles.Padding = 'none'; 
tiles.TileSpacing = 'none'; 

for i = 1:length(files)
    nexttile
    imshow(files{i});
end
%%
close all force; 

files = {   'top_20_pct_stim_consonants_name_left_hemi_left_view.png',
            'top_20_pct_stim_consonants_name_right_hemi_right_view.png',
            'top_20_pct_stim_vowel_name_left_hemi_left_view.png',
            'top_20_pct_stim_vowel_name_right_hemi_right_view.png',
            'top_20_pct_stim_word_name_left_hemi_left_view.png',
            'top_20_pct_stim_word_name_right_hemi_right_view.png',
            'top_20_pct_onset_consonants_name_left_hemi_left_view.png',
            'top_20_pct_onset_consonants_name_right_hemi_right_view.png',
            'top_20_pct_onset_vowel_name_left_hemi_left_view.png',
            'top_20_pct_onset_vowel_name_right_hemi_right_view.png',
            'top_20_pct_onset_word_name_left_hemi_left_view.png',
            'top_20_pct_onset_word_name_right_hemi_right_view.png'};
left_idx = 1:2:length(files);         
left_files = files(left_idx); 

figure();
tiles = tiledlayout(2,3);
tiles.Padding = 'none'; 
tiles.TileSpacing = 'none'; 

for i = 1:length(left_files)
    nexttile
    imshow(left_files{i});
end

right_idx = 2:2:length(files);         
right_files = files(right_idx); 
figure();
tiles = tiledlayout(2,3);
tiles.Padding = 'none'; 
tiles.TileSpacing = 'none'; 

for i = 1:length(right_files)
    nexttile
    imshow(right_files{i});
end        









%%
imshow('print01.jpg')
% imshow('print01.jpg')

%% 
figure
imshow('test.png')