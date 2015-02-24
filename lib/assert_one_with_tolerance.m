function [ t ] = assert_one_with_tolerance( x, tol )
if ~exist('tol','var')
  tol = 10*eps
end
assert( x < 1 + tol );
assert( x > 1 - tol );
end
