function [ t ] = assert_zero_with_tolerance( x, tol )
if ~exist('tol','var')
  tol = 10*eps;
end
assert( abs(x) < tol, 'tolerance exceeded: %e', abs(x) );
end
