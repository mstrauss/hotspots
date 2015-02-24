function [ X,r ] = rmove( obj, X, rmax, distribution )
% RMOVE performs a random move of point X (edgelist format) for random
% amount uniformly between interval [0,rmax) along a random direction.

if ~exist('distribution','var')
  distribution = 'uniform';
end

switch distribution
  case 'uniform'
    r = 2*rmax*rand(size(X,1), 1)-rmax;  % [-rmax, rmax)
  case 'normal'
    r = normrnd(0, rmax, size(X,1), 1);
  otherwise
    error('unknown distribution: %s', distribution);
end

% find random increments
X(:,2) = X(:,2) + r;  % position on edge

% we moved out of the edge
out = X(:,2) < 0 | X(:,2) > obj.L(X(:,1));
if any(out)
  X(out,:) = move_fix( obj, X(out,:) );
end
end
