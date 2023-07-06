classdef classLabelsClass
    properties
        sub_num
        sub_num_formal
        
        word_name
        word_name_cat
        word_name_ohe
        
        consonants_name
        consonants_name_cat
        consonants_name_ohe
        
        vowel_name
        vowel_name_cat
        vowel_name_ohe
    end
    properties(SetAccess = private)
        filelist
    end
    
    methods
        % Constructor
        function obj = classLabelsClass(sub_num)
            obj.sub_num = sub_num; 
            obj.sub_num_formal = strcat('S', string(sub_num));
            persistent filelist 
            if ~isfile('filelist.mat')
                filelist = match_trial_to_behavior();
            else
                filelist = load('filelist.mat', 'filelist').filelist; 
            end
                
            obj.filelist = filelist; 
            sub_trial_data = filelist(filelist.subject == sub_num, 'trials');
            sub_trial_data = sub_trial_data.trials{1}; 
            obj.word_name = table(sub_trial_data.word_name, 'VariableNames', {'word_name'}); 
            obj.word_name_cat = table(categorical(sub_trial_data.word_name), 'VariableNames', {'word_name'}); 
            obj.word_name_ohe = onehotencode(obj.word_name_cat);

            obj.consonants_name = table(sub_trial_data.consonants_name, 'VariableNames', {'consonants_name'}); 
            obj.consonants_name_cat = table(categorical(sub_trial_data.consonants_name), 'VariableNames', {'consonants_name'}); 
            obj.consonants_name_ohe = onehotencode(obj.consonants_name_cat);
            
            obj.vowel_name = table(sub_trial_data.vowel_name, 'VariableNames', {'vowel_name'}); 
            obj.vowel_name_cat = table(categorical(sub_trial_data.vowel_name), 'VariableNames', {'vowel_name'}); 
            obj.vowel_name_ohe = onehotencode(obj.vowel_name_cat);
        end
        
    end
end
            
    
    