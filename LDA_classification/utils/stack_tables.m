function stacked = stack_tables(varargin)
    % STACK_TABLES  Stack tables, inserting NaN values where columns differ.
    %   STACKED = STACK_TABLES(T1) returns T1.
    %   STACKED = STACK_TABLES(T1, T2) joins T2 with T1. If VariableNames (column names) differ between them, the different columns are added and filled with NaN.
    %   STACKED = STACK_TABLES(T1, T2, T3, ... TN) joins T2...TN with T1, substituting NaN where necessary.

    stacked = table;

    for arg_idx = 1:nargin
        t_temp = varargin{arg_idx};
        stacked_col_missing = setdiff(t_temp.Properties.VariableNames, stacked.Properties.VariableNames);
        t_temp_col_missing = setdiff(stacked.Properties.VariableNames, t_temp.Properties.VariableNames);
        stacked = [stacked array2table(nan(height(stacked), numel(stacked_col_missing)), 'VariableNames', stacked_col_missing)];
        t_temp = [t_temp array2table(nan(height(t_temp), numel(t_temp_col_missing)), 'VariableNames', t_temp_col_missing)];
        stacked = [stacked; t_temp];
    end

end
