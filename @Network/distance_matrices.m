function [ B ] = distance_matrices( obj )
% Calculcate weighted (by distance) adjacency matrix.
% B: weighted (euklidian distance between vertices) adjaceny matrix

nn = numnodes(obj.graph);
[s,t] = findedge(obj.graph);
B = sparse( s, t, obj.graph.Edges.Weight, nn, nn );

end
