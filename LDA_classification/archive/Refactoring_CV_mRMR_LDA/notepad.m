%% ToDo:
% refactor ovr_auc_compiled

for class_label = p.class_names
    class_label = class_label{1}
all_class_sub_cv_mRMR_LDA_data.(class_label)
all_class_sub_cv_mRMR_LDA_OvR_data.(class_label)
all_class_sub_indiv_elect_cv_mRMR_LDA_data.(class_label)
all_class_sub_indiv_elect_cv_mRMR_LDA_OvR_data.(class_label)
end

