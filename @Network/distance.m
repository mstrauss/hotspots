function d = distance( obj, p1, p2 )
% DISTANCE calculates the distance along the shortest path between given
% points along the network.
% p1 = [e1, d1]: edge index & distance along this edge for point one
% p2 = [e2, d2; ...]: same for point two but can be MATRIX of points

%warning('Obsolete function');

e1 = p1(:,1);  e2 = p2(:,1);
d1 = p1(:,2);  d2 = p2(:,2);

% length of edge 1
L1 = obj.L(e1);

% same for edge 2
L2 = obj.L(e2);

% % more assertions
% assert( d1 >= 0 );
% assert( all( d2 >= 0 ) );
% assert( d1 <= L1 );
% assert( all( d2 <= L2 ) );

% find vertex indices
[i1,j1] = obj.graph.findedge(e1);
[i2,j2] = obj.graph.findedge(e2);

% minimum distance between all combinations of vertices:
d = zeros( size(e2) );
s = e1==e2;
d(s) = abs( d2(s) - d1 );
i2 = i2(~s);  j2 = j2(~s);
d2 = d2(~s);  L2 = L2(~s);
if nnz(~s)
  d(~s) = min([ ...
    obj.D(i1, i2)' + d1 + d2, ...      % from start to start
    obj.D(i1, j2)' + d1 + L2 - d2, ... % from start to end
    obj.D(j1, i2)' + L1 - d1 + d2, ... % from end to start
    obj.D(j1, j2)' + L1 - d1 + L2 - d2 ], [], 2);
end
end
