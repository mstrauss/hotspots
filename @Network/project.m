function [ C, D, E ] = project( obj, coords )
%PROJECT Projects coords onto network given be edge start vectors
%edge_x and edge length vectors edge_d.
% C: cartesian coordinate represenation of projected coords
% D: displacement (row) vectors
% E: edge/network representation of projected coords, like follows:
%    first column has the edge id,
%    second column the (euclidian) distance of target point from start
%      vertex (of given edge)
%
% NB: Projection is done onto network EDGES, i.e. a projection can happen
% only onto points which also have an edge attached.

% number of input points
N = size(coords, 1);

% number of coordinate dimensions
nd = size(coords, 2);

% number of edges
M = obj.numedges;

% edge vectors
[XD,X1,~,V1,V2] = obj.XD;

% edge lengths
L = obj.L;

% vertex coordinates
XY = obj.vertex_coordinates;

assert( M == size(XD,1) );
assert( nd == size(X1,2) );
assert( nd == size(XD,2) );

% projected output coordinates C and displacement vectors D
C = zeros(N,nd);
D = zeros(N,nd);
E = zeros(N,2);

% norm of a matrix of row vectors
norm = @(X) sqrt(sum(X.^2,2));

% dot product of matrices of row vectors
vdot = @(a, b) sum( bsxfun( @times, a, b ), 2 );

% M-by-nd matrix ; normalized edge_d vectors
[d_hat, norms_d] = normalize(XD);

% calculate closest poinst
[closest_points, closest_points_d] = obj.nearest_vertex( coords );

for i=1:N
  s = coords(i,:);
  
  % length of vector b
  len_b = vdot( d_hat, s ) - vdot( d_hat, X1 );
  b = bsxfun( @times, d_hat, len_b );
  
  delta = bsxfun( @plus, X1, b );
  delta = bsxfun( @minus, delta, s );
  delta_norms = norm(delta);
  [~,J] = sort(delta_norms);
  % indices of valid projections
  % validity depends on the length of b: valid <=> b >= 0 AND b <= |d|
  valid = len_b >= 0 & len_b <= norms_d;
  valid = find( valid(J) );
  if isempty(valid)
    warning('projection failed for index %d, choosing nearest vertex.', i);
    D(i,:) = XY( closest_points(i), : ) - s;
    E(i,:) = map_point(i);
  else
    valid = valid(1);
    valid = paren(paren(1:M,J),valid);
    D(i,:) = delta(valid,:);
    if closest_points_d(i) < delta_norms(valid)
      % disp(['Found a closer vertex for point ', num2str(i)]);
      try
        cp = closest_points(i);
        [E(i,1), E(i,2)] = map_point(i);
        D(i,:) = XY( cp, : ) - s;
      catch
        warning(['Tried to map mark %d [%f x %f] to network point %d [%f x %f], but failed.\n',...
          'Using usual edge mapping.'], ...
          i, s(1), s(2), cp, XY(cp,1), XY(cp,2) );
        [E(i,1), E(i,2)] = map_edge(i);
      end
    else
      [E(i,1), E(i,2)] = map_edge(i);
    end
  end
  C(i,:) = s + D(i,:);
  
  
end

  function [edge_index, d] = map_edge( i )
    edge_index = valid;
    x1 = X1(valid,:);
    x2 = XY(V2(valid),:);
    d = dot( s + D(i,:) -x1 , (x2-x1)/norm(x2-x1) );
  end

  function [edge_index, d] = map_point( i )
    point_index = closest_points(i);
    edge_index = find( V1 == point_index );
    if isempty(edge_index)
      edge_index = find( V2 == point_index );
      if isempty(edge_index)
        error('this point has no edge and is invalid');
      else
        edge_index = edge_index(1);
        d = L(edge_index);
      end
    else
      edge_index = edge_index(1);
      d = 0;
    end
  end

end
