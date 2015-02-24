function [ Z ] = shepard_interpolation( obj, X, t, power )
% Calculate the SHEPARD INTERPOLATION of DATA usind a power law with given
% power as radial basis function.
% X: the coordinates
% t: the targets, i.e. data values t(X)
% Ref: Numerical Recipes, ch. 3.7.3, Shepard Interpolation.

% validate input
if ~exist('power','var'); power = 2; end
[n,d] = size(X);
assert( d == 2 );
assert( all( numel(t) == n ) );
assert( all( size(power) == 1 ) );

% for each data point, calculate the RBF on the grid
Z = zeros( obj.Nx, obj.Ny );
p = -power/2;
for j = 1:obj.Nx*obj.Ny
  [xi,yi] = ind2sub( size(obj), j );
  x = obj.xc(xi);
  y = obj.yc(yi);
  phi = 0;
  yphi = 0;
  for i = 1:n
    % matrix of radial distances
    r2 = ( X(i,1)-x ).^2 + ( X(i,2)-y ).^2;
    if r2 == 0
      Z(j) = t(i);
      break
    else
      % our radial basis function (RBF)
      rbf = r2 .^ p;
      phi = phi + rbf;
      yphi = yphi + t(i) * rbf;
    end
  end
  if Z(j) == 0; Z(j) = yphi / phi; end
end
end
