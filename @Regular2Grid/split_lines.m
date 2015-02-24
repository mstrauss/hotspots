function [E, GI, GJ, DATA, OI, A, V ] = split_lines( obj, network )
% SPLIT_LINES along the grid.
%   All lines are split along the grid lines so that each edge is
%   completely contained inside a grid cell.
% INPUT:
%   network: original network which is to be splitted
% OUTPUT:
%   E: new lines in edgelist format [x y dx dy]
%   GI, GJ: grid cell indices - (1,1) at bottom left
%   DATA: the corresponding metadata
%   OI: the original (unsplitted) edge index
%   A: the updated adjacency matrix
%   V: with correspoding vertex coordinates

% fix input
metadata = network.edgedata;
lines = network.edges;

% validate input
if ~exist('metadata','var')
  metadata = struct;
else
  if ~istable(metadata)
    error('metadata must be a table')
  end
  DATA = table;   % output edgedata
  for name = metadata.Properties.VariableNames
    n = name{1};
    if length(metadata.(n)) ~= size(lines,1)
      error('metadata members must have length equal to no. of lines')
    end
  end
end

% sort lines
lines = normalize_lines( lines );
% find intersections with grid
X = obj.intersectLineSegments( lines );
% translate to grid indices
[xi, xj] = obj.coord2sub( X(:,1), X(:,2) );
XI = [xi,xj];
% init. output variables (dyn. arrays, as we do not know the final amount
% of elements yet)
E = [];  GI = [];  GJ = [];  OI = [];

M = [];  % multiplicity for expanding adjacency matrix
V = zeros(0,2);  % new vertex coordinates
AI = [];  AJ = [];  % indices of new adjaceny matrix

for lineIndex = unique(X(:,3))'  % lines spanning multiple boxes
  line = lines(lineIndex,:);
  line_x0 = line(1:2);  % start point
  line_x1 = line_x0 + line(3:4);  % end point
  % intersection point indices with grid of this line
  I = X(:,3)==lineIndex;
  %yLine = Y(I,:);
  % number of intersections for this line
  n = nnz(I);
  % the relevant (for this line) intersection points
  % (should be sorted!)
  xi = XI(I,:); x = X(I,1:2);
  % locate this edge in adjacency matrix & update multiplicy of adj. matrix
  % expansion
  adjI = network.V1(lineIndex);
  adjJ = network.V2(lineIndex);
  M = [M; adjI adjJ n];
  VI = add_vertex( [line_x0; x; line_x1] );  % VI: vertex indicies
  % for each line segment (i.e. for each box),
  % we add the splitted line segment to the output vectors
  assert( size(xi,1) == n );
  for ls = 1:(n+1)
    % find the box index and length of this line segment
    if ls==1  % one point at the beginning
      [i,j] = obj.coord2sub( line_x0(:,1), line_x0(:,2) );
      E = [E; line_x0, (x(1,:)-line_x0)];
      AI = [AI; VI(1)];  AJ = [AJ; VI(2)];
    elseif ls==n+1  % one point at the end
      [i,j] = obj.coord2sub( line_x1(:,1), line_x1(:,2) );
      E = [E; line_x1, (x(end,:)-line_x1)];
      AI = [AI; VI(end)];  AJ = [AJ; VI(end-1)];
    else  % two points somewhere in-between
      y = xi(ls-1:ls,:);
      E = [E; x(ls-1,1:2), (diff( x(ls-1:ls,1:2) ))];
      AI = [AI; VI(ls)];  AJ = [AJ; VI(ls+1)];
      m = round(mean(y));
      i = m(1); j = m(2);
    end
    GI = [GI; i];
    GJ = [GJ; j];
    OI = [OI; lineIndex];
  end
end

% these lines are contained completely in a box
in_a_box = ~ismember(1:size(lines,1), unique(X(:,3)));
for lineIndex = find(in_a_box)
  line = lines(lineIndex,:);
  line_x0 = line(1:2);  % start point
  VI = add_vertex( [line_x0; line_x0 + line(3:4)] );  % VI: vertex indicies
  AI = [AI; VI(1)];  AJ = [AJ; VI(2)];
  [i,j] = obj.coord2sub( line_x0(:,1), line_x0(:,2) );
  E = [E; line];
  GI = [GI; i];
  GJ = [GJ; j];
  OI = [OI; lineIndex];
end

% resort results
[E, sortindex] = normalize_lines(E);
GI = GI(sortindex);
GJ = GJ(sortindex);
OI = OI(sortindex);

% supply metadata
for name = metadata.Properties.VariableNames
  n = name{1};
  DATA.(n) = metadata.(n)(OI,:);
end

% build new adjacency matrix
n = length(OI);
A = sparse( AI, AJ, ones(length(AI),1), n,n);

function Locb = add_vertex( coords )
  [Lia, Locb] = ismember( coords, V, 'rows' );
  Vlen = size(V,1);
  V = [V; coords(~Lia,:)];
  Locb(~Locb) = Vlen + (1:nnz(~Lia));
end
end
