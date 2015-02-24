function [ B ] = calc_incidence_matrix( obj )
% CALC_INCIDENCE_MATRIX calculates the incidence matrix of the network.

% assert( ~obj.directed );  % (this impl.) only works for undirected networks

B = sparse( obj.N, obj.M );  % rows: vertices, columns: edges

% non-zero elements
[r,c]=find(obj.graph.adjacency);

% find corresponding edges
for i = 1:length(r)
  edge = find( obj.graph.Edges.EndNodes(:,1) == r(i) & obj.graph.Edges.EndNodes(:,2) == c(i) );
  assert( ~isempty(edge) );
  assert( length(edge)==1 );
  B( r(i), edge ) = 1;
  B( c(i), edge ) = 1;
end
