classdef PCA < Domain.Domain
  %% Domain.PCA provides component variables from base variables
  
  properties (Constant)
    VARPREFIX = 'comp';
  end
  
  properties
    % pca inputs
    base_variables
    weight_variable    % an optional variable that is used as weight
    weights            % the evaluated weight variable
    numpoints          % the number of evaluation points
    evaledgecoords     % the edge coordinates (incl network) of the evaluation points
    pca_options        % additional PCA options
    % pca results
    coeff
    score
    latent
    tsquared
    explained
    mu
  end
  
  methods
    function obj = PCA( variables, network, numpoints, weight_variable, pca_options )
      obj.base_variables = AbstractModel.cell2struct( variables );
      assert( isa(network,'Network') );
      obj.numpoints = numpoints;
      [evaldata, obj.evaledgecoords] = obj.variables_to_array( network );

      if exist('weight_variable','var') && ~isempty(weight_variable)
        if ~isa(weight_variable,'Variable')
          error('weights must be a Variable object');
        end
        obj.weight_variable = weight_variable;
      else
        obj.weight_variable = [];
      end
      obj.weights = variable_to_weights( obj, obj.evaledgecoords );
      
      if exist('pca_options','var')
        obj.pca_options = pca_options;
      else
        obj.pca_options = {};
      end
      
      obj.perform_pca( evaldata )
    end
    
    function data = gridded( obj, variable, grid )
      %% calculate the scores of the pca at given grid
      component = strsplit( variable.name, obj.VARPREFIX );
      component = str2double(component(2));
      data = obj.gridded@Domain.Domain(variable, grid);
      nVars = numel(obj.base_variables);
      for i = 1:nVars
        var = obj.base_variables(i);
        data = data + obj.coeff(i,component) * ...
          ( var.gridded( grid ) - obj.mu(i) );
      end
    end
    
    function data = edge_value( obj, variable, edge_coordinates )
      %% calculate the scores of the pca at given edge_coordinates
      component = strsplit( variable.name, obj.VARPREFIX );
      component = str2double(component(2));
      obj.edge_value@Domain.Domain(variable, edge_coordinates);
      data = zeros(edge_coordinates.numpoints, 1);
      nVars = numel(obj.base_variables);
      for i = 1:nVars
        var = obj.base_variables(i);
        data = data + obj.coeff(i,component) * ...
          ( var.edge_value( edge_coordinates ) - obj.mu(i) );
      end
    end
    
    function data = value( obj, variable, X, Y )
      error('implemenation missing');
    end
    
    function bb = bounding_box( obj, variable )
      bb = [NaN NaN; NaN NaN];
    end
    
    function tab = coefftab( obj )
      colNames = obj.coeffnames;
      tab = array2table( obj.coeff, 'RowNames', {obj.base_variables.name}, 'VariableNames', colNames );
    end
    
    function c = coeffnames( obj )
      c = arrayfun( @(i) sprintf('%s%02d',obj.VARPREFIX,i), 1:numel(obj.base_variables), 'UniformOutput', false );
    end
    
    function d = evaldata( obj )
      d = obj.score * obj.coeff' + repmat(obj.mu, size(obj.score,1), 1);
    end
    
    function plotbasevars(obj)
      opts = {'title','PCA Base Variables'};
      if obj.weights
        [N,edges]=histcounts(obj.weights,3);
        groups = discretize( obj.weights, edges );
        opts = [opts {'group',groups}];
      end
      plot_scatter(obj.evaldata,{obj.base_variables.name},false,[],opts{:})
    end
    
    function plotpcavars(obj)
      opts = {'title','PCA Coefficiens vs. Base Variables'};
      if obj.weights
        [N,edges]=histcounts(obj.weights,3);
        groups = discretize( obj.weights, edges );
        opts = [opts {'group',groups,'Y',obj.evaldata,'Ylabels',{obj.base_variables.name}}];
      end
      plot_scatter(obj.score,obj.coeffnames,false,[],opts{:})
    end
  end
  
  methods (Access=private)
    function [array, points] = variables_to_array( obj, net )
      points = EdgeCoordinates( net, net.rpois(obj.numpoints) );
      array = [];
      for i = 1:numel(obj.base_variables)
        var = obj.base_variables(i);
        if iscategorical(var)
          warning('Skipping categorical variable "%s"', var.name);
          continue;
        end
        col = var.edge_value( points );
        array = [array col];
      end
    end
    
    function weights = variable_to_weights( obj, points )
      if isempty(obj.weight_variable)
        weights = false;
        return
      end
      var = obj.weight_variable;
      if iscategorical(var)
        error('Categorical variable "%s" not supported.', var.name);
      end
      weights = var.edge_value( points );
    end
    
    function perform_pca( obj, array )
      opts = obj.pca_options;
      if ~isempty(obj.weight_variable)
        opts = [opts {'Weights', obj.weights}];
      end
      [obj.coeff, obj.score, obj.latent, obj.tsquared, ...
        obj.explained, obj.mu] = pca(array, opts{:});
    end
  end
end
