function [ edgelist, marked_network ] = rpois( obj, n )
% RPOIS Generate a homogenous poisson field with n points.

% build cumulative sum of edge lengths
sL = cumsum(obj.L);

% total network length
netL = sum(obj.L);

% generate random numbers
r = rand(n, 1) * netL;

% find index of edge
I = arrayfun( @(x) paren( find( sL > x ), 1 ), r );

% find positions on that edge
p = sL(I) - r;
% make sure that all positions are within edge length
assert( all( p < obj.L(I) ) );

% contruct edge list
edgelist = [I p];

if nargout > 1
  % construct marked network
  marked_network = MarkedNetwork( obj, edgelist );
end
end
