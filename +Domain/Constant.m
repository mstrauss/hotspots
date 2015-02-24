classdef Constant < Domain.Domain
  
  properties
  end
  
  methods
    function data = gridded( obj, variable, grid )
      obj.gridded@Domain.Domain(variable, grid);
      data = ones(size(grid));
    end
    
    function data = edge_value( obj, variable, edge_coordinates )
      obj.edge_value@Domain.Domain(variable, edge_coordinates);
      data = ones(edge_coordinates.numpoints, 1);
    end
    
    function data = value( obj, variable, X, Y )
      data = ones(numel(X), 1);
    end
    
    function bb = bounding_box( obj, variable )
      bb = [NaN NaN; NaN NaN];
    end
  end
  
end
