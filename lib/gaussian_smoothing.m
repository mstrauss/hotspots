function [ K, N, Nx, Ny, gx, gy ] = gaussian_smoothing( x0, y0, grid_width, kernel_function, bw, Z, gx, gy )
%% GAUSSIAN_SMOOTHING
% x0, y0: coordinates of data
% Z: values of data
%%%
% Notes: DENS_EPS is faster than DENS_DIRECT

% number of data points
N = length(x0);
assert( N == length(y0) );

if (~exist('Z', 'var') || any(isnan(Z(:))) )
  warning('Z contains NaN. Using ones instead');
  Z = ones(1, N, 'double');
end

assert( all(size(Z) == [1 N]), ...
  'Z must be a row vector with length of input data' );

if (~exist('gx', 'var'))
  lx = round(len(x0)*0.05);
  gx = min(x0)-lx:grid_width:max(x0)+lx;
end

% number of query points
Nx = length(gx);

if (~exist('gy', 'var'))
  ly = round(len(y0)*0.05);
  gy = min(y0)-ly:grid_width:max(y0)+ly;
end
Ny = length(gy);

bw2 = 2*bw^2;    % two times bandwidth squared
switch kernel_function
  case 'exp'
    z = bw;
    kernel = @(x,y) exp( -sqrt(x.^2+y.^2)/bw );
  case 'gauss'
    z = sqrt(2*pi)*bw;
    kernel = @(x,y) exp( -(x.^2 + y.^2)/bw2 );
  case 'box'
    kernel = @(x,y) heaviside( bw2 - x.^2 - y.^2 );
  otherwise
    error(['invalid kernel function: ', kernel_function]);
end

[gxx,gyy]=ndgrid( gx, gy ); % grid, auf dem der Kernel berechnet wird
K = zeros(Nx,Ny);
tmp = 0;
with_cutoff = false;
if with_cutoff
  for k=1:N
    cutoff = bw*Z(k);
    x = gxx-x0(k); y = gyy-y0(k);
    filter = abs(x)<cutoff & abs(y)<cutoff;
    kern = kernel( x(filter), y(filter) );
    %tmp = max( tmp, min(kern) );
    K(filter) = K(filter) + Z(k) * kern;
    %K = K + Z(k) * kernel(x,y);
  end
else
  for k=1:N
    x = gxx-x0(k); y = gyy-y0(k);
    K = K + Z(k) * kernel(x,y);
  end
end
% das Ergebnis ist NICHT normalisiert und es sollte sum(K) == sum(Z)
K = K' / z * grid_width^2;
end
