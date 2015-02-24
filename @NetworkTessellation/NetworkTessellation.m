classdef NetworkTessellation < handle
  %NETWORKTESSELLATION Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    weights
    edge_coordinates  % obj of class EdgeCoordinates
    s, t
  end
  
  methods
    function obj = NetworkTessellation( numcells, network, edge_distance_tuples )
      obj.weights = zeros(numcells,1);
      if exist('network','var') && exist('edge_distance_tuples','var')
        obj.edge_coordinates = EdgeCoordinates( network, edge_distance_tuples );
      end
    end
    
    function W = sum( obj )
      % return the sum of weights
      W = sum(obj.weights);
    end
    
    function cell = findcell( obj, edge_indices, edge_distances )
      assert( obj.validate_edge_coordinates( edge_indices, edge_distances ), 'invalid edge coordinates' );
      % FIXME
      error('implementation missing');
    end
  end
  
end
