classdef (Abstract) AbstractModelFit < handle
  
  properties
    parent         % there may be a parent fit where this
    % is a simulation result of
  end
  
  properties (SetAccess = immutable)
    name           % a nice readable name of this fit
    %method         % method string
    model          % a ref. to the parent model
    theta          % mean fit parameters
    meta_names     % names of meta variables
    results        % struct array of fitting results
    resultsMatrix  % as nice matrix
    %resultsTable   % for pretty printing
    configvars     % fitting options
    variates       % a struct array of input data used for fitting
    % (e.g. covariates at data & quadrature points)
    realizations   % realizations of model in edgelist format
    % (i.e. edge id & distance from start vertex)
  end
  
  methods
    function obj = AbstractModelFit( name, model, results, configvars, realizations, variates )
      %% constructor
      obj.name = name;
      %obj.method = method;
      assert( isa(model,'AbstractModel') );
      assert( isstruct(results) );
      assert( isstruct(configvars) );
      assert( isstruct(variates) );
      obj.model = model;
      obj.results = results;
      obj.configvars = configvars;
      obj.variates = variates;
      obj.realizations = realizations;
      
      % collect results
      obj.resultsMatrix = zeros( length(results), ...
        numel(results(1).theta) + numel(results(1).meta ) );
      for i = 1:length(results)
        r = results(i);
        m = numel(r.theta);
        obj.resultsMatrix(i,1:m) = r.theta;
        obj.resultsMatrix(i,m+1:end) = r.meta;
      end
      obj.meta_names = r.meta_names;
      obj.theta = mean(obj.resultsMatrix(:,1:end-obj.number_of_meta_parameters),1);
    end
    
    function v = varnames( obj )
      %% returns the variable names of the fitted model
      v = obj.results(1).fitglm_output.CoefficientNames;
    end
    
    function m = number_of_meta_parameters( obj )
      %% return the number of meta paramters
      if iscell( obj.meta_names )
        m = numel(obj.meta_names);
      else
        m = 1;
      end
    end
    
    function tab = resultsTable( obj, H0_zero_vector )
      if ~exist('H0_zero_vector','var')
        H0_zero_vector = false;
      end
      
      % results table
      if size(obj.resultsMatrix,1) > 1
        if isempty( obj.parent ) || H0_zero_vector
          [z,p] = zScore(obj.resultsMatrix);
        else
          assert( isrow(obj.parent.theta) );
          [z,p] = zScore(obj.resultsMatrix, [obj.parent.theta zeros(1,obj.number_of_meta_parameters)]);
        end
      else
        z = NaN(size(obj.resultsMatrix));
        p = z;
      end
      tab = table( ...
        mean(obj.resultsMatrix,1)', ...
        std(obj.resultsMatrix,[],1)', ...
        z', p', ...
        'VariableNames', {'Mean'; 'SD'; 'Z'; 'p_value'}, 'RowNames', ...
        [obj.varnames obj.meta_names]);
    end
    
  end
end
