classdef GridCollection < handle
  %GridCollection manages a collection of Regular2Grids used for
  %antialiasing.
  
  properties (SetAccess = immutable)
    G        % the combined grid
    N        % number of included grids
    Grids    % cell array of included grids
    bfx, bfy % blow factors, i.e. with what factors do you need to multiply
             % the combined grid's ex, ey to get subgrid's ex, ey.
  end
  
  methods (Static)
    Grids = displaced_grids( base_grid, Ngrids )
    [seq, bfx, bfy] = offset_sequence( Ngrids )
  end
  
  methods
    function obj = GridCollection( base_grid, Ngrids )
      assert(isa(base_grid,'Regular2Grid'));
      obj.Grids = GridCollection.displaced_grids( base_grid, Ngrids );
      obj.N = numel(obj.Grids);
      if obj.N > 1
        [obj.G, obj.bfx, obj.bfy] = obj.combine();
      else
        obj.G = base_grid;
        obj.bfx = 1;
        obj.bfy = 1;
      end
    end
    
    function bb = bounding_box( obj )
      bb = obj.G.bounding_box;
    end
    
    function bbi = bounding_box_ind( obj )
      bbi = obj.G.bounding_box_ind;
    end
  end
  
  methods (Access = private)
    [resulting_grid, bfx, bfy] = combine( obj )
    [resulting_data] = prepare( obj, data, bbi )
  end
  
end
