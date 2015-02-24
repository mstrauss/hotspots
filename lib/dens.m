function [ K, N, Ny, gx, gy ] = dens( x0, y0, Nx, bw, Z, gx, gy )
%DENS Summary of this function goes here
%   Detailed explanation goes here

N = length(x0);
assert( N == length(y0) );

if (~exist('Z', 'var') | isnan(Z) )
  Z = ones(1, N, 'double');
end

if (~exist('gx', 'var'))
  gx = linspace(min(x0),max(x0),Nx);
else
  Nx = length(gx);
end

if (~exist('gy', 'var'))
  Ny = round(Nx/len(x0)*len(y0));
  gy = linspace(min(y0),max(y0),Ny);
else
  Ny = length(gy);
end

gauss2d = @(x,y,sx,sy) exp(-x.^2 / sx^2/2) .* exp(-y.^2 / sy^2/2) / 2/pi/sqrt(sx^2*sy^2);
gauss2dbw = @(x,y) gauss2d(x,y,bw,bw);

[gxx,gyy]=ndgrid( gx, gy ); % grid, auf dem der Kernel berechnet wird
K = zeros(Nx,Ny);
for k=1:N
  K = K + Z(k) * gauss2dbw( gxx-x0(k), gyy-y0(k) );
end
K = K'; %/N;

end
