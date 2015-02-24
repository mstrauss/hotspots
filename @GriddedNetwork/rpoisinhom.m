function [ edgelist, accept_rate, marked_network ] = rpoisinhom( obj, n, dist_fun )
% RPOISINHOM Generate a inhomogenous poisson field with n points.
%   dist_fun is a the 2d proability density (need not be normalized) given
%   on dist_grid (Regular2Grid object) from which the random field is
%   sampled using the naïve rejection method.

%% validate input
assert( all( size( dist_fun ) == size(obj.grid) ) );
net = obj.network;

% build cumulative sum of edge lengths
sL = cumsum(net.L);

% total network length
netL = sum(net.L);

% find maximum value of sampling distribution
fMax = max( dist_fun(:) );

% number of random positions to generate per sweep
m = 10*n;

% number of found points
N = 0;

% found points
II = zeros(n,1); pp = zeros(n,1);

% statistics to calculate acceptance rate
acc = 0; total = 0;

while N < n
  % generate random positions on the network
  rx = rand(m, 1) * netL;

  % find index of edge
  I = arrayfun( @(x) paren( find( sL > x ), 1 ), rx );

  % find positions on that edge
  p = sL(I) - rx;

  % find xy coordinates
  rxy = net.xy_from_edge( I, p );

  % find ij indices on grid
  [ri, rj] = obj.grid.coord2sub( rxy(:,1), rxy(:,2) );

  % generate random numbers to decide acceptance/rejection
  rf = rand(m, 1) * fMax * 1.01;
  % times 1.01, so that we have also some probability to hit the maximum

  % keep only positions rx where rf is < sampling_distribution
  fxy = dist_fun( sub2ind( size(dist_fun), ri, rj ) );
  sel = rf < fxy;
  new = nnz(sel);
  acc = acc + new;  total = total + m;
  if new > n-N
    new = n-N;
    sel = paren( find(sel), 1:new );
  end
  II(1+N:N+new) = I(sel);
  pp(1+N:N+new) = p(sel);
  N = N + new;
end

% contruct edge list
edgelist = [II pp];

if nargout > 1
  % calculate acceptance rate
  accept_rate = acc/total;
end

if nargout > 2
  % construct marked network
  marked_network = MarkedNetwork( obj, edgelist );
end

end
