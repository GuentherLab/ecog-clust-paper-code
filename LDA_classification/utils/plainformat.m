function s = plainformat(s)
s = regexprep(s, '[\\\^\_]','\\$0');
end