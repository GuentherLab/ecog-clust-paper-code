function out = plaintext(str_or_cell)

if iscell(str_or_cell)
    str = str_or_cell{1};
else
    str = str_or_cell;
end

c = strsplit(str, '_');
str = strjoin(c);
    
out = regexprep(str, '(^|\s).', '${upper($0)}');
end

    