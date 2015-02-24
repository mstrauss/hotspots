function [z, p] = zScore( vector, H0_mu )
% calculate z-score and p-value for given vector under the null hypothesis
% that the vector mean is mu (double sided).
if ~exist('H0_mu','var')
  H0_mu = zeros(1,size(vector,2));
end
assert( all( size(H0_mu) == [1 size(vector,2)] ) );
if ~isvector(vector)
  [z,p] = arrayfun( @(i) zScore(vector(:,i), H0_mu(i)), 1:size(vector,2) );
else
  [mu, sigma] = normfit(vector);
  z = (mu-H0_mu)/sigma;
  p = 2*normcdf( abs(z), 0, 1, 'upper' );
end
end
