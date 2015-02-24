classdef Model < AbstractModel
  %Network.MODEL describes a model living on a network.
  
  properties (SetAccess = immutable)
    network    % here lives the model
  end  
 
  methods
    
    function obj = Model( network, type, vars )
      obj@AbstractModel( type, vars );
      assert( isa(network,'Network') );
      obj.network = network;
    end
    
  end

  methods (Access=private)
  end
end
