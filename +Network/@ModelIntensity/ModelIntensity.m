classdef ModelIntensity < NetworkTessellation
  %MODELINTENSITY Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (SetAccess=immutable)
    prob
  end
  
  methods
    function obj = ModelIntensity( edge_coordinates, probabilities )
      assert( isa( edge_coordinates, 'EdgeCoordinates' ) );
      obj@NetworkTessellation( edge_coordinates.net, edge_coordinates.edge_network );
      assert( numel(probabilities) == edge_coordinates.numpoints );
      assert( sum(probabilities)-1 < 1e-10, 'probabilites do not sum up to 1' );
      obj.prob = probabilities;
    end
    
    function N = numcells( obj )
      N = obj.numpoints;
    end
    
    function values = value( obj, edge_indices, edge_distances )
      % find a probability on given edge at given distance
      values = obj.prob( obj.findcell( edge_indices, edge_distances ) );
    end
    
  end
  
end
