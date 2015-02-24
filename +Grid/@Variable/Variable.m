classdef Variable < handle
  
  properties (SetAccess=immutable)
    data      % the gridded variable data
    varname   % the gridded variable name
    variable  % handle to the original variable
  end
  
  methods
    function obj = Variable( varname, griddata, variable )
      assert( isa(variable,'Variable') );
      obj.data = griddata;
      obj.varname = varname;
      obj.variable = variable;
    end
  end
  
end
