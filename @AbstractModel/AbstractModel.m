classdef (Abstract) AbstractModel < handle
  %AbstractModel describes an abstract model (i.e. without domain).
  
  properties (SetAccess = immutable)
    type       % the type of the model
    variables  % cell array of variables (including intercept)
  end
  
  methods (Static)
    function types = types()
      % returns the valid types
      types = {'loglinear'};
    end
    
    function s = cell2struct( vars )
      s = Variable.empty(numel(vars),0);
      for i = 1:numel(vars)
        v = vars{i};
        if isa(v,'Variable')
          s(i) = v;
        else
          s(i) = Variable(v);
        end
      end
    end
  end

  methods (Abstract)
    Vq = edge_value( obj, variable, edge_coordinates );
  end
  
  methods    
    function obj = AbstractModel( type, vars )
      assert( ismember(type, AbstractModel.types) );
      obj.type = type;
      obj.variables = AbstractModel.cell2struct( vars );
      obj.validate_variables
    end
    
    function varnames = varnames( obj, idx )
      varnames = {obj.variables.name};
      if exist('idx','var'); varnames = varnames(idx); end
    end
    
    function cat = categorical( obj, idx )
      cat = [obj.variables.categorical];
      if exist('idx','var'); cat = cat(idx); end
    end
    
  end

  methods (Access=private)
    function validate_variables( obj )
      if strcmp(obj.type,'loglinear')
        % obj.validate_log_transform
      end
    end
    
    function validate_log_transform( obj )
      for var = obj.variables
        if ( ~var.categorical && ~isa(var.domain,'Domain.Constant') ) || ...
            ~isa(var.transform,'function_handle') || ~strcmp('log',func2str(var.transform))
          warning('Variable %s has no logarithmic transform in a log-linear model.\nThis m', var.name );
        end
      end
    end
  end
end
