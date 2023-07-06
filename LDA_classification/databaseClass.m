classdef databaseClass
    properties
        electrodes_database
        subject_class_labels
    end
    methods
        function db = databaseClass(electrodes_database)
            db.electrodes_database = electrodes_database; 
            
            subject_class_labels = table();
            for sub_num = unique(electrodes_database.subject)
                subject_class_labels_temp = table({string(sub_num)}, sub_num, {classLabelsClass(sub_num)},...
                    'VariableNames', {'sub_num_formal', 'sub_num', 'class_labels'});
                subject_class_labels = [subject_class_labels; subject_class_labels_temp]; 
                
            end
           
            db.subject_class_labels = subject_class_labels;
        end
    end
end