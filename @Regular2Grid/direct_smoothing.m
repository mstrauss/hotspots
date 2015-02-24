function [ Vq, N ] = direct_smoothing( obj, X, Y, V, bandwidth, kernel_function )

if ~exist('kernel_function','var'); kernel_function = 'gaussian'; end

[ Vq, N ] = gaussian_smoothing( X, Y, obj.grid_width, kernel_function, ...
  bandwidth, V, obj.xc, obj.yc );
Vq = Vq';
