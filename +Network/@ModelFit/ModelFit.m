classdef ModelFit < AbstractModelFit
  %% Network.ModelFit
  
  properties (SetAccess=immutable)
    marked_network % the fitting objective
  end
  
  methods
    function obj = ModelFit( marked_network, name, model, results, configvars, realizations, variates )
      obj@AbstractModelFit( name, model, results, configvars, realizations, variates );
      assert( isa(marked_network,'MarkedNetwork') );
      obj.marked_network = marked_network;
    end
    
  end
  
end
