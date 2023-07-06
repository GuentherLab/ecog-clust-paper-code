% feat_check
if ~isvarname('asit')
    asit = load('/usr3/graduate/lcj126/Shortcut_ECoG_Preprocessed_LJ/MLCompiledData/FeatureSets/e2000_w50_s25/indiv_elect_feature_set_table.mat').all_subs_indiv_elect_feat_table;
end

truth = [];
for idx = 1:height(edb)
    sub_num = edb.subject(idx);
    edb_type = edb.type(idx);
    edb_electrode = edb.electrode(idx);
    edb_feat = edb.feat_set(idx);
    
    if edb_type == 2
        asit_feat = asit.word_name.stim_feat{asit.sub_num == sub_num};
    elseif edb_type == 1
        asit_feat = asit.word_name.onset_feat{asit.sub_num == sub_num};
    end

    t = asit_feat.(2)(asit_feat.(1) == edb_electrode);
    t = t{1};
    t = t.(1); 
    t = t(1,:);
    
    x = edb_feat{1};
    x = x(1,:);
    
    truth(idx) = isequal(t,x);
    

    
    
end
    
if sum(truth) == height(edb)
    disp('all values verified');
else
    disp('values between tables not equivalent');
end
    
    
    

