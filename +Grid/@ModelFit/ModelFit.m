classdef ModelFit < AbstractModelFit
  %Grid.MODELFIT
  
  properties
  end
  
  methods
    function obj = ModelFit( name, model, results, configvars, realizations, variates )
      obj@AbstractModelFit( name, model, results, configvars, realizations, variates );
    end
  end
  
end
