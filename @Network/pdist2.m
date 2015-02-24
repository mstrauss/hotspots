function [ D ] = pdist2( obj, X, Y )
%% PDIST2 calculates distances between points X and Y along the network.
%% X, Y: edgelist format.

% number of edges
M = size(X,1);
N = size(Y,1);
D = zeros(M,N);
for i = 1:M
  j = 1:N;
  D(i,j) = obj.distance( X(i,:), Y(j,:) );
end
end
