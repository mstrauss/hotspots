function [ tess, unaccounted ] = weights( obj )
% Calculates weights as network length per mark.
% It does so, by iterating over all network edges and assigning to each
% mark accordingly.
% W..vector of weights
% unaccounted...total unaccounted network length
%
% NB: This function implements KIND-OF Voronoi Network tesselation.
% Differences: THIS METHOD WORKS WITH NON-SIMPLE NETWORK LAYOUTS, i.e.
%   WE ALLOW FOR CRITICAL NODES.
% See 2012-Okabe/Sugihara-Spatial Analysis along Networks Statistical and
% Computational Methods, chap. 4.3.

% debugging?
dbg = false;

if dbg
  obj.plot('debug',1);
end

unaccounted = 0;

% fetch the incidence matrix of the network
I = obj.net.graph.incidence;

% no of network edges
N = obj.net.numedges;

% initialize weights, indexed by ordering of sorted_edge_network
% (this is fixed at the end of the function)
%W = zeros( M, 1 );

% sorting the mark coordinates
[sorted_edge_network, original_mark_index] = sortrows(obj.edge_network);

% retrieve network metrics
V1 = obj.net.V1;
V2 = obj.net.V2;
L = obj.net.L;

% resulting tessellation, tuples of [edge_index edge_distance] representing
% network cell boundaries
tess = NetworkTessellation( obj.nummarks );
% we can discern btw. 2 types of network cells: shared and non-shared;  the
% shared cells belong to more than one mark;  the resulting assignments of
% marks to network cells are stored as tuples of [cell

if dbg; old_W_sum = 0; end

% loop over all edges
for edge = 1:N
  if dbg; fprintf('Calculating for edge %d (L=%f):\n', edge, obj.net.L(edge) ); end

  % n1: neighboring marks of start vertex
  n1 = I( V1(edge), : );
  n1(edge) = 0;
  mark_distances_1 = zeros(nnz(n1),1);
  mark_ids_1 = zeros(nnz(n1),1);
  n1i = find(n1);
  for j = 1:nnz(n1)
    visited_edges = false(1,N);
    visited_edges(edge) = true;
    [mark_distances_1(j), mark_ids_1(j)] = mark_distance( n1i(j), V1(edge) );
    %mark_distances_1(j) = obj.graph.Nodes.D(V1(i));
    %mark_ids_1(j) = obj.graph.Nodes.T(V1(i));
  end
  % n2: neighboring marks of end vertex
  n2 = I( V2(edge), : );
  n2(edge) = 0;
  mark_distances_2 = zeros(nnz(n2),1);
  mark_ids_2 = zeros(nnz(n2),1);
  n2i = find(n2);
  for j = 1:nnz(n2)
    visited_edges = false(1,N);
    visited_edges(edge) = true;
    [mark_distances_2(j), mark_ids_2(j)] = mark_distance( n2i(j), V2(edge) );
    %mark_distances_2(j) = obj.graph.Nodes.D(V2(i));
    %mark_ids_2(j) = obj.graph.Nodes.T(V2(i));
  end

  if dbg
    fprintf( 'mark_distances_1: %f\n', mark_distances_1 );
    fprintf( 'mark_ids_1: %d\n', mark_ids_1 );
    fprintf( 'mark_distances_2: %f\n', mark_distances_2 );
    fprintf( 'mark_ids_2: %d\n', mark_ids_2 );
    fprintf( '----\n' );
  end
  
  % marks on currently enumerated edge
  marks_on_edge_index = sorted_edge_network(:,1) == edge;
  no_of_marks_on_edge = nnz(marks_on_edge_index);
  mark_distances = sorted_edge_network(marks_on_edge_index,2);
  assert( all( mark_distances < L(edge) ) );
  marks_in_edge_network = find(marks_on_edge_index);
  
  if no_of_marks_on_edge > 0
    % assign start part weights to marks
    assignable_length = mark_distances(1);
    if assignable_length > 0
      sel = mark_distances_1 < assignable_length;
      assign_weights( tess, mark_ids_1(sel), mark_distances_1(sel), false )
    end
    % assign end part weights to marks
    assignable_length = L(edge) - mark_distances(end);
    if assignable_length > 0
      sel = mark_distances_2 < assignable_length;
      assign_weights( tess, mark_ids_2(sel), mark_distances_2(sel), true )
    end
    for j = 1:no_of_marks_on_edge-1
      assignable_length = mark_distances(j+1)-mark_distances(j);
      if assignable_length > 0
        local_weights = 0.5 * ones(2,1);
        local_mark_ids = [marks_in_edge_network(j); marks_in_edge_network(j+1)];
        %assert(~isempty(local_mark_ids));
        %W(local_mark_ids) = W(local_mark_ids) + assignable_length * local_weights/sum(local_weights);
        tess.weights(local_mark_ids) = tess.weights(local_mark_ids) + assignable_length * local_weights/sum(local_weights);
      end
    end
  else  % no_of_marks_on_edge == 0
    assignable_length = L(edge);
    local_weights = [ 1./(assignable_length+mark_distances_1); 1./(assignable_length+mark_distances_2) ];
    valid_weights = local_weights>0;
    local_weights = local_weights(valid_weights);
    if isempty(local_weights)
      unaccounted = unaccounted + assignable_length;
    else
      local_mark_ids = [ mark_ids_1; mark_ids_2 ];
      local_mark_ids = local_mark_ids(valid_weights);
      % account for possible duplicate mark ids
      %assert(~isempty(local_mark_ids));
      [s,si] = sort(local_mark_ids);
      [a,~,c] = unique(s);
      asum = assignable_length * local_weights(si)/sum(local_weights);
      %W(a) = W(a) + accumarray( c, asum );
      tess.weights(a) = tess.weights(a) + accumarray( c, asum );
    end
  end
  
  % checkup
  if dbg
    assert( abs( sum(tess.weights) + unaccounted - old_W_sum - obj.net.L(edge) ) < 100000*eps );
    old_W_sum = sum(tess.weights) + unaccounted;
  end
end

% check that for the total network length is accounted for
assert( abs( (sum(tess.weights) + unaccounted - sum(obj.net.L))/sum(obj.net.L) ) < 1e-10 );

% fix ordering of weights
tess.weights = tess.weights(original_mark_index);
if unaccounted > 0
  warning('Disconnected network? %f unaccounted for.', unaccounted);
end

  function [D, Marks, L] = mark_distance_vectorized( edge_index_vector, start_vertex )
    D = inf( length(edge_index_vector), 1 );
    Marks = zeros( length(edge_index_vector), 1 );
    L = zeros( length(edge_index_vector), 1 );
    for k = 1:length(edge_index_vector)
      [D(k), Marks(k), L(k)] = mark_distance( edge_index_vector(k), start_vertex);
    end
  end

  function [d, mark, path_len] = mark_distance( edge_index, start_vertex )
    %%MARK_DISTANCE returns the distance d to the closest mark, starting
    %%from given edge and vertex.  Also returns the mark index and the
    %%total path length. BUT: does not search on the starting edge!
    %assert( isscalar(start_vertex) );
    visited_edges(edge_index) = true;
    % marks on that edge
    marks_on_edge_index = sorted_edge_network(:,1) == edge_index;
    table_offset = find( marks_on_edge_index, 1 ) - 1;
    no_of_marks_on_edge = nnz(marks_on_edge_index);
    path_len = L( edge_index );
    if no_of_marks_on_edge > 0
      marks_on_sorted_edge_network = sorted_edge_network(marks_on_edge_index,:);
      if start_vertex == V1(edge_index)
        [d, idx] = min( marks_on_sorted_edge_network(:,2) );
      else
        %assert( start_vertex == V2(edge_index) );
        [tmp, idx] = max( marks_on_sorted_edge_network(:,2) );
        d = L(edge_index) - tmp;
      end
      mark = table_offset + idx;
    else
      % need to look further
      if start_vertex == V1(edge_index)
        new_start_vertex = V2(edge_index);
      else
        %assert( start_vertex == V2(edge_index) );
        new_start_vertex = V1(edge_index);
      end
      neighboring_edges = I( new_start_vertex, : );
      neighboring_edges(visited_edges) = 0;
      neighboring_edges_index = find(neighboring_edges);
      [mark_dists, mark_ids, path_lens] = mark_distance_vectorized( neighboring_edges_index, new_start_vertex );
      [d, minIndex] = min([Inf; mark_dists]);
      if minIndex > 1;  mark = mark_ids(minIndex-1); else mark = 0; end
      path_len = path_len + sum(path_lens);
      % update inf to the subnetwork length
      d = d + L( edge_index );
    end
  end

  function assign_weights( tess, mark_ids, mark_distances, from_end )
    % assignable_length lies either on the beginning or the end
    % (iff from_end=true) of the given edge
    local_weights = [ 1/assignable_length; 1./(assignable_length+mark_distances) ];
    valid_weights = local_weights>0;
    local_weights = local_weights(valid_weights);
    if from_end
      local_mark_ids = marks_in_edge_network(end);
    else
      local_mark_ids = marks_in_edge_network(1);
    end
    local_mark_ids = [ local_mark_ids; mark_ids ];
    local_mark_ids = local_mark_ids(valid_weights);
    % account for possible duplicate mark ids
    %assert(~isempty(local_mark_ids));
    [s,si] = sort(local_mark_ids);
    [a,~,c] = unique(s);
    asum = assignable_length * local_weights(si)/sum(local_weights);
    %W(a) = W(a) + accumarray( c, asum );
    tess.weights(a) = tess.weights(a) + accumarray( c, asum );
  end
end


