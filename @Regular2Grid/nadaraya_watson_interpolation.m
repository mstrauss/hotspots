function [ Vq, N ] = nadaraya_watson_interpolation( obj, X, Y, V, bandwidth )

if ~exist('kernel_function','var'); kernel_function = 'gaussian'; end

%%% slow:
[Xq, Yq] = obj.ndgrid;
[ Vq, N ] = nadaraya_watson_interpolation( X(:), Y(:), V(:), Xq(:), Yq(:), bandwidth );
Vq = reshape( Vq, obj.Nx, obj.Ny )';
