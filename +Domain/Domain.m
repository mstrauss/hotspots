classdef (Abstract) Domain < handle
  %DOMAIN represents a variable domain, i.e. where a variable is defined.
  %  This class defines the interface for the subclasses of Domain.
  
  properties
  end
  
  methods
    function data = gridded( obj, variable, grid )
      obj.assert_variable( variable );
      assert( isa(grid,'Regular2Grid') );
      if nargout > 0; data = zeros(size(grid)); end
    end
    
    function Vq = edge_value( obj, variable, edge_coordinates )
      obj.assert_variable( variable );
      assert( isa(edge_coordinates,'EdgeCoordinates'), 'object of class EdgeCoordinates req.' );
      if nargout > 0
        [Xq,Yq] = edge_coordinates.XY;
        Vq = obj.value( variable, Xq, Yq );
      end
    end
    
    function Vq = value( obj, variable, Xq, Yq )
      obj.assert_variable( variable );
      assert( isnumeric(Xq) );
      assert( isnumeric(Yq) );
      assert( numel(Xq) == numel(Yq) );
      if nargout > 0; Vq = zeros(numel(Xq), 1); end
    end
  end
  
  methods (Abstract)
    bb = bounding_box( obj, variable );
  end
  
  methods (Access=protected)
    function assert_variable( obj, variable )
      assert( isa(variable,'Variable') );
      assert( variable.domain == obj );
    end
  end
  
end
