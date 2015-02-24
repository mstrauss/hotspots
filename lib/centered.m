function [ X ] = centered( X, X0 )
if ~exist('X0','var')
  X0 = X;
end
mX0 = mean(X0,1);
assert( all( size(mX0) == [1, size(X,2)] ) );
% calculate bin-wise mean
X = X ./ repmat(mX0,size(X,1),1);
end

