function [ S, K ] = convolution_smoothing( obj, sigma, data )
% CONVOLUTION_SMOOTHING smoothes DATA by convoluting them with a Gaussian
% kernel with covariance matrix SIGMA.  DATA must have the same size as the
% grid. SIGMA must be a symmetric 2x2-matrix.

%% Clear NaN and Inf
data(isnan(data))=0;
data(isinf(data))=0;

%% Input Validation
if ~exist('sigma','var'); sigma = [1 0; 0 1]; end
assert( all( size(data) == [obj.Nx obj.Ny] ) );
assert( all( size(sigma) == [2 2] ) );

sigma_det = det(sigma);

%% Construct kernel meshgrid, depending on bandwith
expansion_factor = 20;
mesh_width = ceil( expansion_factor * sqrt(norm(sigma*[1 0]')/obj.ex) );
mesh_height = ceil( expansion_factor * sqrt(norm(sigma*[0 1]')/obj.ey) );
mesh_x = -mesh_width/2*obj.ex:obj.ex:mesh_width/2*obj.ex;
mesh_y = -mesh_height/2*obj.ey:obj.ey:mesh_height/2*obj.ey;
[X,Y] = meshgrid( mesh_x, mesh_y );
assert( all( size(X) == [mesh_height+1, mesh_width+1] ) );
assert( all( size(Y) == [mesh_height+1, mesh_width+1] ) );

%% Caclulate kernel on meshgrid
Kfun = @(x,y) exp( -([x y]/sigma*[x;y])/2 );
K = arrayfun( Kfun, X, Y ) * ( obj.ex * obj.ey / (2*pi*sqrt(sigma_det)));
%assert(all(isodd(size(K))));
Keps = abs(sum(K(:))-1);
if Keps > 1e-5
  warning('Keps=%.2g > 1e-5: Kernel seems to have bad quality.', Keps)
end
K = K ./ sum(K(:));

%% Perform convolution
S = conv2( data, K, 'same' );
end
