function [ Vq ] = parzen_estimation( obj, X, bandwidth )
%% PARZEN_ESTIMATION
% applies kernel density estimation aka Parzen estimation to data points X
% on the grid.
%%%
% See 2006-Bishop-Pattern Recognition and Machine Learning, p. 123.

assert( size(X,2) == 2, 'need 2-dimensional input (row) vectors in X');

% number of data points
N = size(X,1);

Vq = reshape( ...
  ksdensity( X, obj.coordinates, 'Bandwidth', bandwidth), ...
  obj.Nx, obj.Ny );

end

