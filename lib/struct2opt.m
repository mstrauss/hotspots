function c = struct2opt(s)
% source: http://stackoverflow.com/a/15013309
fname = fieldnames(s);
fval = struct2cell(s);
c = [fname, fval]';
c = c(:);
