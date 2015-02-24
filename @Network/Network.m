classdef Network < handle
  %NETWORK class
  
  properties (SetAccess = immutable)
    % calculated properties
    B       % weighted (distances) upper adjacency matrix
    D       % matrix of shortest distances between vertices
    Dfun    % what kind of distance function has been used
    %I       % incidence matrix, rows: vertices, cols: edges
  end
  
  properties (SetAccess = private)
    graph   % the graph object
  end
  
  methods
    function obj = Network(vertex_coordinates, options)
      if ~isstruct(options)
        error('need an options struct argument');
      end
      
      EdgeTable = table();
      NodeTable = table();
      
      % handle node coordinates
      NodeTable.XY = vertex_coordinates;
      NodeTable.X = vertex_coordinates(:,1);
      NodeTable.Y = vertex_coordinates(:,2);
      
      if isfield(options,'edgelist')
        V1 = [options.edgelist.v1];
        V2 = [options.edgelist.v2];
      elseif isfield(options,'adjacency')
        assert(issymmetric(options.adjacency), 'A needs to be symmetric');
        [V1,V2] = find(triu(options.adjacency));
      elseif isfield(options,'weights')
        W = options.weights;
        if ~iscolumn(W); W=W'; end
      else
        error('need either "edgelist" or "adjacency" option');
      end
      
      if ~iscolumn(V1); V1=V1'; end
      if ~iscolumn(V2); V2=V2'; end
      
      X1 = vertex_coordinates(V1,:);
      X2 = vertex_coordinates(V2,:);
      XD = X2 - X1;
      assert( ~any(all( XD==0, 2 )), sprintf('edge %d has zero length', find( all( XD==0, 2 ),1 )) );
      EdgeTable.EndNodes = [V1 V2];
      if exist('W','var')
        EdgeTable.Weight = W;
      else
        EdgeTable.Weight = vnorm(XD);
      end
      
      if isfield(options,'edgedata')
        error( 'edgedata support is broken' );
        if istable(options.edgedata)
          for f = options.edgedata.Properties.VariableNames
            f = f{:};
            if ismember(f, {'EndNodes','Weight'})
              %  assert(all(all( options.edgedata.(f) == EdgeTable.(f) )));
            else
              EdgeTable.(f) = options.edgedata.(f);
            end
          end
        else
          error('illegal edgedata type');
        end
      end
      
      [EdgeTable, iToUniq] = unique(EdgeTable);
      assert( max(hist(iToUniq,numel(EdgeTable)))==1, 'Oops, hit duplicate edges.' );
      
      % build graph
      obj.graph = graph( EdgeTable, NodeTable );
      
      % validate graph
      assert( all( obj.graph.Edges.Weight == sqrt( ...
        diff(obj.graph.Nodes.X(obj.graph.Edges.EndNodes),[],2).^2 + ...
        diff(obj.graph.Nodes.Y(obj.graph.Edges.EndNodes),[],2).^2 ) ) );
      
      if isfield(options,'distance')
        obj.Dfun = options.distance;
      else
        obj.Dfun = 'dijkstra';
      end
      
      switch obj.Dfun
        case 'dijkstra'
          obj.D = obj.graph.distances;
        case 'euclidean'
          if obj.N > 1e4
            error('distance not calculated. too much data.');
          else
            obj.D = pdist(obj.vertices);
          end
        otherwise
          error('distance must be specified');
      end
      
      obj.B = obj.distance_matrices;
      %obj.I = obj.calc_incidence_matrix;
    end
    
    function x = x(obj, varargin)
      x = obj.graph.Nodes.X;
      x = x(varargin{:});
    end
    
    function y = y(obj, varargin)
      y = obj.graph.Nodes.Y;
      y = y(varargin{:});
    end
    
    function N = N(obj)
      warning('obsolete function');
      % return number of vertices
      N = obj.graph.numnodes;
    end
    
    function N = numnodes(obj)
      N = obj.graph.numnodes;
    end
    
    function M = M(obj)
      warning('obsolete function');
      % return number of edges
      M = obj.graph.numedges;
    end
    
    function M = numedges(obj)
      M = obj.graph.numedges;
    end
    
    function gp = plot(obj, plot_labels, varargin)
      if exist('plot_labels','var') && plot_labels
        nodeLabels = arrayfun(@(i) {sprintf('V.%d', i)}, (1:obj.numnodes)');
        edgeLabels = arrayfun(@(i) {sprintf('E.%d', i)}, (1:obj.numedges)');
        plot(obj.graph,'XData',obj.graph.Nodes.X,'YData',obj.graph.Nodes.Y,...
          'NodeLabel', nodeLabels, 'EdgeLabel', edgeLabels, varargin{:})
        axis xy equal
      else
        gp = plot(obj.graph,'XData',obj.graph.Nodes.X,'YData',obj.graph.Nodes.Y,...
          'NodeLabel', [], 'Marker', 'none', varargin{:})
      end
    end
    
    function d = diameter( obj )
      % diameter of largest component
      d = max(obj.D(obj.D<Inf));
    end
    
    function [I, L] = get_edge( obj, i, j )
      % find edge index and edge length from given vertex indices
      % index of edge
      I = obj.graph.findedge(i,j);
      if I == 0  % not in graph
        I = NaN;
        L = 0;
      else
        % length of edge
        L = obj.graph.Edges.Weight(I);
        % assert correct direction
        assert( obj.graph.Edges.EndNodes(I,1) == i, 'invalid direction' );
        assert( obj.graph.Edges.EndNodes(I,2) == j );
      end
    end
    
    function X = xy_from_edge( obj, edge_indices, edge_distances )
      % xy_from_edge maps edge coordinates [index, offset] to xy coordinates.
      assert( obj.validate_edge_coordinates( edge_indices, edge_distances ), 'invalid edge coordinates' );
      
      % calculate coordinates
      [v1, v2] = obj.vertices_for_edges( edge_indices );
      x1 = obj.vertex_coordinates(v1);
      x2 = obj.vertex_coordinates(v2);
      X = normalize( x2-x1 ) .* repmat(edge_distances', 2, 1)' + x1;
    end
    
    function filtered_network = edge_filter(obj, edge_indices)
      if islogical(edge_indices)
        assert( numel(edge_indices) == obj.M );
      else
        assert( all(edge_indices > 0) && all(edge_indices <= obj.M) );
      end
      filtered_network = Network( obj.vertex_coordinates, struct(...
        'distance', obj.Dfun, ...
        'edgelist', obj.edgelist( edge_indices ) ) );
    end
    
    function [v1,v2] = vertices_for_edges( obj, edge_indices )
      % return vertex indices for given edge indices
      if islogical(edge_indices); edge_indices = find(edge_indices); end
      [v1,v2] = obj.graph.findedge(edge_indices);
    end
    
    function edgelist = edgelist( obj, edge_indices )
      if ~exist('edge_indices','var')
        [v1,v2] = obj.graph.findedge;
      else
        [v1,v2] = vertices_for_edges( obj, edge_indices );
      end
      edgelist = struct('v1', v1', 'v2', v2');
    end
    
    function V1 = V1( obj, idx )
      % start vertices of all edges
      if ~exist('idx','var')
        V1 = obj.graph.Edges.EndNodes(:,1);
      else
        V1 = obj.graph.Edges.EndNodes(idx,1);
      end
    end
    
    function V2 = V2( obj, idx )
      % terminal vertices of all edges
      if ~exist('idx','var')
        V2 = obj.graph.Edges.EndNodes(:,2);
      else
        V2 = obj.graph.Edges.EndNodes(idx,2);
      end
    end
    
    function XY = vertex_coordinates( obj, vertex_indices )
      if ~exist('vertex_indices','var'); vertex_indices = 1:obj.numnodes; end
      XY = [obj.graph.Nodes.X(vertex_indices) obj.graph.Nodes.Y(vertex_indices)];
    end
    
    function L = L(obj, edge_indices)
      % return length of given edges
      if ~exist('edge_indices','var'); edge_indices = 1:obj.numedges; end
      L = obj.graph.Edges.Weight(edge_indices);
    end
    
    function length = length(obj)
      length = sum(obj.graph.Edges.Weight(:));
    end
    
    function [xd,x1,x2,v1,v2] = XD(obj, edge_indices)
      % displacement vectors; X2 = X1 + XD
      if ~exist('edge_indices','var'); edge_indices = 1:obj.numedges; end
      [v1, v2] = obj.vertices_for_edges( edge_indices );
      x1 = obj.vertex_coordinates(v1);
      x2 = obj.vertex_coordinates(v2);
      xd = x2-x1;
    end
    
    function edges = edges(obj, edge_indices)
      if ~exist('edge_indices','var'); edge_indices = 1:obj.numedges; end
      [xd,x1] = obj.XD(edge_indices);
      edges = [x1,xd];
    end
    
    function edgedata = edgedata(obj)
      edgedata = obj.graph.Edges;
    end
    
    function bb = bounding_box( obj )
      bb = [min(obj.graph.Nodes.X) min(obj.graph.Nodes.Y);
        max(obj.graph.Nodes.X) max(obj.graph.Nodes.Y)];
    end
    
    function tf = validate_edge_coordinates( net, E, D )
      % validate edge_indices (E) and edge distances (D)
      tf = all( E > 0 & E <= net.numedges ) && all( D >=0 & D <= net.L(E) );
    end
    
    function tf = validate_nodes( net, nodeidx )
      tf = all( nodeidx > 0 & nodeidx <= net.numnodes );
    end
    
    function net = split_edges( obj, splitpoints )
      % return a new network with split edges
      E = splitpoints(:,1);
      D = splitpoints(:,2);
      assert( obj.validate_edge_coordinates( E, D ), 'invalid split points' );
      edgelist = obj.edgelist;
      newnode = obj.numnodes + 1;
      [exv1,exv2] = obj.vertices_for_edges( splitpoints(:,1) );
      % keep track of remapped nodes and edges
      nodemap = 1:obj.numnodes;
      for i = 1:numel(E)
        v1 = nodemap( exv1(i) );
        v2 = nodemap( exv2(i) );
        edge = find( edgelist.v1==v1 & edgelist.v2==v2 );
        % remap
        edgelist.v2( edge ) = newnode;
        nodemap(v1) = newnode;
        % append
        edgelist.v1 = [edgelist.v1 newnode];
        edgelist.v2 = [edgelist.v2 v2];
        % update trackers
        nodemap(v1) = newnode;
        % update counter
        newnode = newnode + 1;
      end
      net = Network( [obj.vertex_coordinates; obj.xy_from_edge(E,D)], ...
        struct('edgelist', edgelist) );
    end
    
    function net = delete_nodes( obj, nodeidx )
      % return a new network with given nodes deleted
      assert( obj.validate_nodes( nodeidx ), 'invalid nodes' );
      
      % track remapped nodes
      nodemap = 1:obj.numnodes;
      % map removed nodes to new edges
      node2edge = obj.graph.incidence ~= 0;
      
      % what nodes to keep
      keepnodes = setdiff( nodemap, nodeidx );
      
      % identify to-be-deleted edges
      endnodes = obj.graph.Edges.EndNodes;
      weights  = obj.graph.Edges.Weight;
      [tf,loc] = ismember( endnodes, nodeidx );
      
      % output edgelist & weights
      edgelist = struct;  edgelist.v1 = [];  edgelist.v2 = [];
      weightsout = [];
      
      % loop over edges
      for edgeidx=1:obj.numedges
        % keep edge
        if ~any(tf(edgeidx,:))
          edgelist.v1 = [edgelist.v1; endnodes(edgeidx,1)];
          edgelist.v2 = [edgelist.v2; endnodes(edgeidx,2)];
          weightsout = [weightsout; weights(edgeidx)];
        else
          A = tf(edgeidx,1);  B = tf(edgeidx,2);
          if xor( A, B )
            if A
              rewire( nodeidx(loc(edgeidx,1)) );
            else
              rewire( nodeidx(loc(edgeidx,2)) );
            end
          else  % rewire both nodes
            rewire( nodeidx(loc(edgeidx,:)) );
          end
        end
      end
      
      % validate
      assert( numel(weightsout) ==  numel(edgelist.v1), 'wrong number of weights' );
      assert( abs( obj.length - sum(weightsout) ) < 1e-10, 'wrong net length' );
      
      % construct new network
      net = Network( obj.vertex_coordinates, ...
        struct('edgelist', edgelist, 'weights', weightsout) );
      
      function rewire( nodeidx )
        nodeidx = nodemap(nodeidx)
        posnodeidx = nodeidx( nodeidx > 0 );
        if ~all(nodeidx); return; end  % nothing to do
        switch length(posnodeidx)
          case 0
            warning('skipping edge %s, (L=%f)',edgeidx,obj.L(edgeidx));
            return
          case 1  % rewire single node
            k = obj.graph.degree(posnodeidx);
            assert( k == 2, 'only k=2 supported' );
            n = obj.graph.neighbors(posnodeidx);
          case 2  % rewire edge
            n = cell2mat( arrayfun( @(idx) obj.graph.neighbors(idx), posnodeidx, 'UniformOutput', false ) );
          otherwise
            error('rewiring multiple nodes not supported');
        end
        n = setdiff( setdiff( n(:), posnodeidx), find(~nodemap) );
        if numel(n) < 2
          % assign weight to first neighboring edge
          %w = weights(edgeidx);
          %assgnedge = find( node2edge(:,edgeidx),1 );
          %weightsout(assgnedge) = weightsout(assgnedge) + w;
        else
          % create new edge
          s = n(1);  t = n(2);
          edgelist.v1 = [edgelist.v1; s];
          edgelist.v2 = [edgelist.v2; t];
          % update trackers
          nodemap(nodeidx)=false;
          oldedges = node2edge(posnodeidx,:);
          node2edge(:,oldedges) = 0;
          [~,oldedges]=find(oldedges);  oldedges=unique(oldedges);
          % update weights
          weightsout = [weightsout; sum(weights(oldedges))];
        end
      end
    end
  end
  
  methods (Static)
    function net = createFromLines( lines )
      assert( isnumeric(lines) );
      assert( size(lines,2) == 4 );
      nLines = size(lines,1);
      % create a Network from lines, i.e. tuples of [x0, y0, x1, y1]
      startX = lines(:,1:2);
      endX = lines(:,3:4);
      [vertex_coordinates,~,loc] = unique( [startX; endX], 'rows' );
      edgelist.v1 = loc(1:nLines);
      edgelist.v2 = loc((1+nLines):end);
      net = Network( vertex_coordinates, ...
        struct('edgelist',edgelist) );
    end
  end
  
end
