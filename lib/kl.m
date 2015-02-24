function [ kl ] = kl( p, q )
%KL calculates Kullback-Leibler divergence for given input densities

if isstruct(p); p = p.z; end
if isstruct(q); q = q.z; end
p = p ./ sum(p(:));
q = q ./ sum(q(:));
assert( all( p(:) >= 0 ) );
assert( all( q(:) >= 0 ) );
kl = p .* log(p./q);
kl(isnan(kl))=0;
if any(isinf(kl(:)))
  %warning('Some values Inf.  Resetting to zero.');
  kl(isinf(kl))=0;
end
end
