classdef EdgeCoordinates < handle
  %EDGECOORDINATES represents a set of edge coordnates,
  % tuples of [edge_index, edge_offset] to identify network locations.
  % It is a lightweight MarkedNetwork.
  
  properties (SetAccess=immutable)
    net          % a network
    edge_network % edgelist representation: [edge_index distance]
  end
  
  methods
    function obj = EdgeCoordinates( net, edge_coordinates )
      assert( isa(net, 'Network'), 'need a Network as first paramter' );
      obj.net = net;

      assert( isnumeric(edge_coordinates) );
      assert( size(edge_coordinates,2) == 2 );
      E = edge_coordinates(:,1);
      D = edge_coordinates(:,2);
      
      assert( obj.validate_edge_coordinates( E, D ), 'invalid edge coordinates' );
      obj.edge_network = edge_coordinates;
    end
    
    function E = E(obj)
      E = obj.edge_network(:,1);
    end
    
    function D = D(obj)
      D = obj.edge_network(:,2);
    end
    
    function M = M(obj)
      warning( 'obsolete function, use "numpoints" instead' );
      % number of points
      M = size(obj.edge_network,1);
    end
    
    function N = numpoints(obj)
      N = size(obj.edge_network,1);
    end
    
    function [X,Y] = XY(obj)
      XY = obj.net.xy_from_edge( obj.E, obj.D );
      X = XY(:,1);
      Y = XY(:,2);
    end
  end
  
  methods (Access = protected)
    function tf = validate_edge_coordinates( obj, edge_indices, edge_distances )
      tf = obj.net.validate_edge_coordinates( edge_indices, edge_distances );
    end
  end
end
