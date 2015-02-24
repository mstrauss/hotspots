function [ Counts, I, J, CountsMatrix ] = boxcount( obj, X, Y )
assert( all( size(X) == size(Y) ) );

% inside grid?
assert( min(X) >= obj.xgMin );
assert( max(X) <= obj.xgMax );
assert( min(Y) >= obj.ygMin );
assert( max(Y) <= obj.ygMax );

CountsMatrix = zeros(obj.Nx, obj.Ny);
[I,J] = obj.coord2sub( X, Y );
for row = 1:numel(I)
  CountsMatrix(I(row), J(row)) = CountsMatrix(I(row), J(row)) + 1;
end

% return list of non-zero counts & coordinates
if ~iscolumn(I); I=I'; end
if ~iscolumn(J); J=J'; end
X = unique([J I], 'rows');
X=sortrows(X);
I = X(:,2);
J = X(:,1);
Counts = CountsMatrix(CountsMatrix>0);

end
