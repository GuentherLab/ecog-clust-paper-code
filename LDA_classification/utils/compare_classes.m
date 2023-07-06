function class_diff_table = compare_classes(c1, c2)
    % COMPARE_CLASSES  Find the difference between two classes.
    %   COMPARE_CLASSES(C1,C2) finds the properties of each class, their values, then
    %   checks for equivalence WRT the other. If they are different values, the
    %   property in question is placed in a table.
    %
    %   Properties unique to C1/C2 are also added, with NaN in place as
    %   the other class' value for that property.
    % formatter ignore 2
    c1_props = properties(c1);
    c2_props = properties(c2);

    c1_unique_props = setdiff(c1_props, c2_props);
    c2_unique_props = setdiff(c2_props, c1_props);

    shared_props = union(c1_props(ismember(c1_props, c2_props)), c2_props(ismember(c2_props, c1_props)));
    class_diff_table = table();

    for shared_prop_idx = 1:length(shared_props)
        shared_prop = char(shared_props(shared_prop_idx));
        c1_currentProperty = getfield(c1, shared_prop);
        c2_currentProperty = getfield(c2, shared_prop);

        if ~isequal(c1_currentProperty, c2_currentProperty)
            class_diff_table = stack_tables(class_diff_table, table({shared_prop}, {c1_currentProperty}, {c2_currentProperty}, ...
                'VariableNames', {'property', 'c1_value', 'c2_value'}));
        end

    end

    for c1_unique_prop_idx = 1:length(c1_unique_props)
        c1_unique_prop = c1_unique_props{c1_unique_prop_idx};
        c1_unique_prop_val = getfield(c1, c1_unique_prop);
        class_diff_table = stack_tables(class_diff_table, table({c1_unique_prop}, {c1_unique_prop_val}, NaN, ...
            'VariableNames', {'property', 'c1_value', 'c2_value'}));
    end

    for c2_unique_prop_idx = 1:length(c2_unique_props)
        c2_unique_prop = c2_unique_props{c2_unique_prop_idx};
        c2_unique_prop_val = getfield(c2, c2_unique_prop);
        class_diff_table = stack_tables(class_diff_table, table({c2_unique_prop}, NaN, {c2_unique_prop_val}, ...
            'VariableNames', {'property', 'c1_value', 'c2_value'}));
    end

    if height(class_diff_table) > 0
        class_diff_table = movevars(class_diff_table, 'property', 'Before', 1);
    end

end
