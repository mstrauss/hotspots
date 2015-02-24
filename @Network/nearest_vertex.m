function [ I, D ] = nearest_vertex( obj, S )
%NEAREST_VERTEX finds the nearest vertex from (row vector) coordinates X
%for each (row) coordinate vector in S.
% Outputs:
%   I: array of indices of closest vertices
%   D: array of distances (to closest vertices)

N = size( S, 1 );
nd = size( S, 2 );
assert( nd == size( obj.vertex_coordinates, 2 ) );

I = zeros( N, 1 );
D = zeros( N, 1 );

for i = 1:N
  distances = sqrt(sum( ( bsxfun( @minus, obj.vertex_coordinates, S(i,:) ) ).^2, 2) );
  closest = find( distances == min(distances) );
  if length(closest) > 1
    % disp(['Warning: ambigious closest point for index ', num2str(i)]);
    closest = closest(1);
  end
  I(i) = closest;
  D(i) = distances( I(i) );
end
end
