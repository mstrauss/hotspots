classdef NetworkEdge < Domain.Domain
  
  properties
    network
  end
  
  methods
    function obj = NetworkEdge( network )
      assert(isa(network, 'Network'));
      obj.network = network;
    end
    
    function data = gridded( obj, variable, grid )
      obj.gridded@Domain.Domain(variable, grid);
      edgedata = obj.network.edgedata.(variable.name);
      
      % apply transform, if applicable
      if isa(variable.transform,'function_handle')
        edgedata = variable.transform(edgedata);
      end
      
      if variable.categorical
        categories = unique(edgedata);
        categories = categories(categories>0);
        data = cell(numel(categories),1);
        for i = 1:numel(categories)
          data{i} = struct();
          level = categories(i);
          data{i}.level = level;
          d = grid.line_density( obj.network.edges(edgedata==level) );
          d(isnan(d))=0;
          data{i}.data = d;
        end
      else
        error('only categorical variables supported at this point');
      end
      
    end
    
    function data = edge_value( obj, variable, edge_coordinates )
      obj.edge_value@Domain.Domain(variable, edge_coordinates);
      assert( edge_coordinates.net == obj.network );
      edgedata = obj.network.edgedata.(variable.name);

      % apply transform, if applicable
      if isa(variable.transform,'function_handle')
        edgedata = variable.transform(edgedata);
      end

      data = edgedata(edge_coordinates.E);
    end
    
    function data = value( obj, variable, X, Y )
      error('to be implemented');
    end
    
    function bb = bounding_box( obj, variable )
      bb = obj.network.bounding_box;
    end
    
  end
  
end
