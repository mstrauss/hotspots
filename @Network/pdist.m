function [ D ] = pdist( obj, X )
%% PDIST calculates distances between points of X along the network.
%% X: edgelist format.

% number of edges
M = size(X,1);
D = zeros(M);
% work through upper triangonal matrix
for i = 1:M
  j = i+1:M;
  D(i,j) = obj.distance( X(i,:), X(j,:) );
end
% assert( all( diag(D) == 0 ) );
D = D + D';
end
