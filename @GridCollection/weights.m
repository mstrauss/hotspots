function [W,K] = weights( obj, subgridindex )
%% this calculates a matrix of Gaussian weights for given subgrid
% the weight is maximal, if subgrid and mainpixel have the same center

G = obj.G;
sG = obj.Grids{subgridindex};

% convolution kernel
m = obj.bfx( subgridindex );
n = obj.bfy( subgridindex );

X = -(m-1)/2:1:m/2;
Y = -(n-1)/2:1:n/2;
YX = combinations(Y,X);
K = reshape( mvnpdf(YX, 0, diag([m n]) ), m, n );

% this is completely the same as in GridCollection#prepare
dx0 = ( ( obj.G.xgMin - obj.Grids{subgridindex}.xgMin )/obj.G.ex );
dy0 = ( ( obj.G.ygMin - obj.Grids{subgridindex}.ygMin )/obj.G.ey );
dx1 = ( ( obj.Grids{subgridindex}.xgMax - obj.G.xgMax )/obj.G.ex );
dy1 = ( ( obj.Grids{subgridindex}.ygMax - obj.G.ygMax )/obj.G.ey );
dx0 = int32(dx0);  dx1 = int32(dx1);
dy0 = int32(dy0);  dy1 = int32(dy1);

W = repmat(K, 1+sG.Nx, 1+sG.Ny );
offsets = obj.offset_sequence();
dx0 = offsets(subgridindex,1); dy0 = offsets(subgridindex,2);
W = W(1+dx0:obj.G.Nx+dx0, 1+dy0:obj.G.Ny+dy0);
assert( all( size(W) == size(G) ) );
end
