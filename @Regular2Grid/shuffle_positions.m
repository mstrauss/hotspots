function [ newS ] = shuffle_positions( obj, S, dist_fun )
%% shuffle positions of struct S within grid
if ~isfield(S,'x') || ~isfield(S,'y')
  error('x,y fields required.')
end

assert( all( size(S.x) == size(S.y) ) );
n = length(S.x);
assert( n > 1 );

if ~exist('dist_fun', 'var')
  dist_fun = ones( obj.Nx, obj.Ny );
end

assert( all( size(dist_fun) == [obj.Nx, obj.Ny] ), 'dist_fun has wrong dimension' )

% find maximum value of sampling distribution
fMax = max( dist_fun(:) );

% number of random positions to generate per sweep
m = 10*n;

% number of found points
N = 0;

% found points
XX = zeros(n,1); YY = zeros(n,1);

% statistics to calculate acceptance rate
acc = 0; total = 0;

while N < n
  % generate random positions
  rx = rand(m, 1) * (obj.xgMax-obj.xgMin) + obj.xgMin;
  ry = rand(m, 1) * (obj.ygMax-obj.ygMin) + obj.ygMin;
  
  % find ij indices on grid
  [ri, rj] = obj.coord2sub( rx, ry );
  
  % generate random numbers to decide acceptance/rejection
  rf = rand(m, 1) * fMax * 1.01;

  % keep only positions rx where rf is < sampling_distribution
  fxy = dist_fun( sub2ind( size(dist_fun), ri, rj ) );
  sel = rf < fxy;
  new = nnz(sel);
  acc = acc + new;  total = total + m;
  if new > n-N
    new = n-N;
    sel = paren( find(sel), 1:new );
  end
  XX(1+N:N+new) = rx(sel);
  YY(1+N:N+new) = ry(sel);
  N = N + new;
end

% return results as row vectors
newS = struct('x', XX', 'y', YY');

end
