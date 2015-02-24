function [ edgelist, accept_rate, marked_network ] = rpoisinhom( model_intensity, n )
% RPOISINHOM Generate a inhomogenous poisson field on the network with n
% points.  model_intensity is a ModelIntensity object from which the random
% field is sampled using the naïve rejection method.

%% validate input
assert( isa( model_intensity, 'ModelIntensity' ) );
net = model_intensity.net;

% build cumulative sum of edge lengths
sL = cumsum(net.L);

% total network length
netL = sum(net.L);

% find maximum value of sampling distribution
fMax = max( model_intensity.prob );

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

  % generate random numbers to decide acceptance/rejection
  rf = rand(m, 1) * fMax * 1.01;
  % times 1.01, so that we have also some probability to hit the maximum

  % keep only positions rx where rf is < sampling_distribution
  fxy = model_intensity.value( I, p );
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
  marked_network = MarkedNetwork( net, edgelist );
end

end
