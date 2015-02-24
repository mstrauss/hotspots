function [ x, y ] = select_curve( X, Y, fun, q )
% takes curves (rows of X and Y) and pointwisely applies fun (e.g. @min or @max)
% and returns the single resulting curve

assert( size(X,2) == size(Y,2) );
assert( size(X,1) == size(Y,1) || size(X,1) == 1 );
n = size(Y,1);
m = size(Y,2);

if ~exist('q','var')
  q = 0.025;
end

if isa(fun, 'function_handle')
  if strcmp(func2str(fun), 'median')
    y = nanmedian(Y,1);
    %I = arrayfun( @(i) find(Y(:,i)==y(i), 1, 'first' ), 1:m );
    I = getindex;
  else
    [y,I] = fun(Y,[],1);
  end
else
  qIndex = ceil( (n+1)*q );
  Ysorted = sort(Y);  % sorted columnwise
  if strcmp(fun, 'lower')
    y = Ysorted(qIndex,:);
  elseif strcmp(fun, 'upper')
    y = Ysorted(end-qIndex+1,:);
  elseif strcmp(fun, 'mean')
    y = nanmean( Y );
  end
  I = arrayfun( @(i) min([ find(Y(:,i)==y(i), 1, 'first') 1]), 1:m );
end
if size(X,1) > 1
  x = arrayfun( @(x,y) X(x,y), I, 1:m );
else
  x = X;
end
[x, I] = sort(x);
y = y(I);

  function tmp = getindex
    for i = 1:m
      tmp = find( Y(:,i) == y(i), 1, 'first' );
      if isempty(tmp)
        tmp = 0;
      end
    end
  end

end
