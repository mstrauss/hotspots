classdef MarkedNetwork < handle
  %MarkedNetwork marks on a network
  
  properties (SetAccess = immutable)
    net          % a network
    edge_network % edgelist representation: [edge_index distance]
    %MD           % matrix of shortest distances between mark
    edge_coordinates_ordering
  end
  
  properties
    filter       % filter marks for given indices
  end
  
  properties (SetAccess = private)
    graph
    distances  % cache for matrix of shortest distances btw. marks
  end

  methods (Static)
    function mn = fromCartesianMarkCoordinates( net, XY )
      [~,~,edges] = net.project(XY);
      mn = MarkedNetwork( net, edges );
    end
  end
  
  methods
    
    function N = N(obj)
      warning('obsolete function');
      % number of marks
      N = nnz(obj.graph.Nodes.Type==2);
    end
    
    function N = nummarks(obj)
      N = nnz(obj.graph.Nodes.Type~=1);
    end
    
    function obj = MarkedNetwork( net, edge_coordinates )
      % MarkedNetwork( net, coords, use_edge_coordinates )
      %   net: a network object
      %   edge_coordinates: coordinates, given as [edge_index edge_offset]
      
      % sort edge coordinates
      [edge_coordinates, edge_coordinates_ordering] = sortrows( edge_coordinates );
      obj.edge_network = edge_coordinates;
      obj.edge_coordinates_ordering = edge_coordinates_ordering;
      
      obj.net = net;
      obj.graph = net.graph;
      obj.graph.Nodes.Type = ones(net.numnodes,1);
      
      obj.add_marks( edge_coordinates );
      
      obj.filter = 1:obj.nummarks;
    end
    
    function coords = coords(obj)
      coords = obj.graph.Nodes.XY(obj.graph.Nodes.Type==2,:);
    end
    
    function x = x(obj, varargin)
      x = obj.graph.Nodes.XY(obj.graph.Nodes.Type==2,1);
      x = x(varargin{:});
    end
    
    function y = y(obj, varargin)
      y = obj.graph.Nodes.XY(obj.graph.Nodes.Type==2,2);
      y = y(varargin{:});
    end
    
    function X = coords_from_edge_network( obj, varargin )
      % calculate mark coordinates from edge_network
      edge_index = obj.edge_network(:,1);
      [xd,x1] = obj.net.XD( edge_index );
      d = normalize( xd );
      X = bsxfun( @times, d, obj.edge_network(:,2) ) + x1;
      X = X(varargin{:});
    end
    
    function gp = plot( obj, varargin )
      % parse options
      argi = 1;
      while argi <= numel(varargin)
        switch varargin{argi}
          case 'plot_labels'
            argi = argi + 1; plot_labels = varargin{argi};
          case 'plot_orig'
            argi = argi + 1; plot_orig = varargin{argi};
          case 'plot_options'
            argi = argi + 1; plot_options = varargin{argi};
          case 'add_marks'
            argi = argi + 1; add_marks = varargin{argi};
          case 'marker_size'
            argi = argi + 1; marker_size = varargin{argi};
          case 'debug'
            argi = argi + 1; debug = varargin{argi};
          otherwise
            error('Invalid option "%s". Syntax changed.', varargin{argi});
        end
        argi = argi + 1;
      end
      % set default options
      if ~exist('plot_labels','var'); plot_labels = false; end
      if ~exist('plot_options','var'); plot_options = '+r'; end
      if ~exist('add_marks','var'); add_marks = false; end
      if ~exist('marker_size','var'); marker_size = 6; end
      if ~exist('line_style','var'); line_style = '-'; end
      if exist('debug','var') && debug
        obj.net.plot(0);
        hold on
        plot_labels = true;
        line_style = 'none';
      end
      
      gp = plot(obj.graph,'XData',obj.graph.Nodes.X,'YData',obj.graph.Nodes.Y);
      gp.Marker = repmat({'none'},1,obj.graph.numnodes);
      marks = find( obj.graph.Nodes.Type == 2 );
      marks = marks( obj.filter );
      gp.Marker( marks ) = {'o'};
      gp.NodeColor = 'red';
      gp.MarkerSize = marker_size;
      gp.LineStyle = line_style;
      if plot_labels
        if debug
          gp.NodeLabel = cellfun( @(c) num2str(c), num2cell(1:obj.graph.numnodes), 'UniformOutput', 0 );
          gp.EdgeLabel = cellfun( @(c) num2str(c), num2cell(1:obj.graph.numedges), 'UniformOutput', 0 );
        else
          gp.NodeLabel = cellfun( @(c) num2str(c), cell(obj.graph.numnodes,1), 'UniformOutput', 0);
          gp.NodeLabel(marks) = cellfun( @(c) num2str(c), num2cell(1:obj.nummarks), 'UniformOutput', 0 );
        end
      else
        gp.NodeLabel = [];
      end
    end
    
    function gp = plot_components( mn )
      gp = plot(mn.graph,'XData',mn.graph.Nodes.X,'YData',mn.graph.Nodes.Y);
      gp.Marker = repmat({'none'},1,mn.graph.numnodes);
      gp.NodeColor = zeros(mn.graph.numnodes,3);
      netnodes = ( mn.graph.Nodes.Type == 1 );
      marknodes = ( mn.graph.Nodes.Type ~= 1 );
      gp.Marker( marknodes ) = {'x'};
      gp.Marker( netnodes ) = {'o'};
      gp.MarkerSize = 6;
      comp = mn.graph.conncomp;
      col = hsv(numel(unique(comp)));
      gp.NodeColor = col( comp', : );
    end
    
    function obj = set.filter(obj,indices)
      obj.filter = indices;
    end
    
    function [n1,n2,s] = nearest_neighbor_distances(obj)
      s = sort(obj.MD,2);
      n1 = s(:,2);
      n2 = s(:,3);
    end
    
    function s = is_simple(obj)
      s = ( size(unique(obj.edge_network,'rows'),1) == obj.nummarks );
    end
    
    function MD = MD(obj)
      if isempty(obj.distances)
        switch obj.net.Dfun
          case 'euclidean'
            obj.distances = squareform(pdist( obj.coords ));
          case 'dijkstra'
            sel = find( obj.graph.Nodes.Type==2 );
            obj.distances = obj.graph.distances(sel,sel);
          otherwise
            error('not implemented');
        end
      end
      MD = obj.distances;
    end
    
    function D = distances_types( obj, t1, t2 )
      % calculate distances between two types of vertices
      sel1 = find( obj.graph.Nodes.Type==t1 );
      sel2 = find( obj.graph.Nodes.Type==t2 );
      D = obj.graph.distances(sel1,sel2);
    end
    
    function add_marks( obj, edge_coordinates, type )
      % add marks of given type to network
      
      % default type
      if ~exist( 'type', 'var' ); type = 2; end
      
      % find coordinates of marks
      edge_index = edge_coordinates(:,1);
      edge_dist = edge_coordinates(:,2);
      xy = obj.net.xy_from_edge( edge_index, edge_dist );
      
      % add marks as new nodes of Type <type>
      N = size(edge_coordinates,1);
      obj.graph = obj.graph.addnode(table(xy, xy(:,1), xy(:,2), repmat(type,N,1), 'VariableNames',{'XY','X','Y','Type'}));
      
      % remember data from old edges
      old_v1 = obj.graph.Edges.EndNodes(edge_coordinates(:,1),1);
      old_v2 = obj.graph.Edges.EndNodes(edge_coordinates(:,1),2);
      old_L = obj.graph.Edges.Weight(edge_coordinates(:,1));
      
      % find ids of new nodes
      marks = ( numnodes(obj.net.graph)+1 : numnodes(obj.net.graph)+N )';
      
      % find the single and the multi-marks per edge cases
      [unique_edges, eci, uei] = unique( edge_coordinates(:,1) );
      % number of marks per unique edge
      nMarksUE = accumarray(uei,1);
      
      % now, handle all single-mark-per-edge cases
      sel = eci(nMarksUE==1);
      if nnz(sel) > 0
        edgeT = repmat( obj.graph.Edges( edge_coordinates(sel,1), : ), 2,1 );
        edgeT.EndNodes = [old_v1(sel) marks(sel); marks(sel) old_v2(sel)];
        edgeT.Weight = [edge_coordinates(sel,2); old_L(sel) - edge_coordinates(sel,2)];
      else
        edgeT = [];
      end
      
      % remaining (unused) mark indices
      marks = setdiff( marks, marks(sel) );
      marksIndex = 1;
      
      % now, handle the multi-mark cases
      sel = nMarksUE > 1;
      unique_edges = unique_edges(sel);  eci = eci(sel);
      for unique_edge_index = 1:numel(unique_edges)
        unique_edge = unique_edges(unique_edge_index);
        old_edge_index = eci(unique_edge_index);
        distances_on_edge = edge_coordinates( edge_coordinates(:,1) == unique_edge, 2 );
        m = size(distances_on_edge,1);
        t = repmat( obj.graph.Edges( unique_edge, : ), m+1,1 );
        Lrem = old_L(old_edge_index);  % how much length remains
        Lgiven = 0;  % how much length given away
        prev_vertex = old_v1( old_edge_index );
        next_vertex = marks(marksIndex);  marksIndex = marksIndex + 1;
        next_seglen = distances_on_edge(1);
        for i = 1:m+1
          % add record to table
          t.EndNodes(i,:) = [prev_vertex next_vertex];
          t.Weight(i) = next_seglen;
          Lgiven = Lgiven + next_seglen;
          prev_vertex = next_vertex;
          if i < m
            next_vertex = marks(marksIndex);  marksIndex = marksIndex + 1;
            next_seglen = distances_on_edge(i+1)-Lgiven;
          elseif i == m
            next_vertex = old_v2(old_edge_index);
            next_seglen = Lrem - Lgiven;
          else
            assert_zero_with_tolerance( Lrem - Lgiven );
          end
        end
        edgeT = [edgeT; t];
      end
      
      obj.graph = obj.graph.addedge(edgeT);
      obj.graph = obj.graph.rmedge( old_v1, old_v2 );
      
      % validation
      assert_zero_with_tolerance( sum(obj.net.graph.Edges.Weight) - sum(obj.graph.Edges.Weight), 1e-10 );
    end
    
    function k = degree( mn, nodeidx, type )
      t1 = mn.graph.Nodes.Type == 1;
      k = arrayfun( @(n) nnz(mn.graph.Nodes.Type( mn.graph.neighbors(n) )==type), nodeidx );
    end
    
    %%% OBSOLETE METHODS
    
    function n = NU(obj)
      warning('MarkedNetwork#NU is obsolete and will be removed');
      % for compatiblity
      n = size(unique(obj.edge_network,'rows'),1);
    end
    
    function coords = expanded_coords(obj)
      warning('MarkedNetwork#expanded_coords is obsolete and will be removed');
      % for compatiblity
      coords = obj.coords;
    end
    
    function en = expanded_edge_network(obj)
      warning('MarkedNetwork#expanded_edge_network is obsolete and will be removed');
      % for compatiblity
      en = obj.edge_network;
    end
    
    function MI = mark_incidence( mn )
      % contract the incidence matrix such that only marks remain
      md = triu(mn.MD);
      g = graph(md+md');
      MI = g.incidence;
    end
    
  end
  
  
  
end
