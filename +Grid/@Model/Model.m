classdef Model < Network.Model
  %%Grid.Model is a model that uses pixelimages to represent the variates.
  
  properties (SetAccess = immutable)
    G          % the grid, where the model is calculated
    Gvariables % the model variables transformed to the grid
  end
  
  methods
    function obj = Model( network, type, vars, grid_width )
      obj@Network.Model( network, type, vars );
      
      % init grid
      [~, obj.G] = load_dataset( 'pop', grid_width );
      
      % init grid data
      for vari = 1:numel(obj.variables)
        var = obj.variables(vari);
        data = var.gridded( obj.G );
        if iscell(data)
          % split categorical variables
          for i = 1:numel(data)
            GX = data{i}.data;
            lvl = data{i}.level;
            if ~ischar(lvl); lvl = num2str(lvl); end
            varname = sprintf('%s_%s', var.name, lvl);
            obj.Gvariables = [obj.Gvariables ...
              Grid.Variable( varname, GX, var )];
          end
        else
          obj.Gvariables = [obj.Gvariables ...
            Grid.Variable( var.name, data, var )];
        end
      end
      
    end
    
    function Vq = edge_value( obj, variable, edge_coordinates )
      %% retrieve variable value at given coordinates in [edgeidx distance] format
      if variable.categorical
        Vq = variable.edge_value( edge_coordinates );
      else
        [X,Y] = edge_coordinates.XY;
        sel = obj.G.coord2ind( X, Y );
        %Vq = interpn( xx, yy, obj.GX(:,:,i), X(:,1), X(:,2), interpolation_method );
        Vq = obj.Gvariables(strcmp(obj.Gvarnames,variable.name)).data(sel);
      end
    end
    
    function grid = grid( obj )
      grid = obj.G;
    end
    
    function names = Gvarnames( obj )
      names = {obj.Gvariables.varname};
    end
    
    function data = GX( obj, varidx )
      data = obj.Gvariables(varidx).data;
    end
        
    function plotvar( obj, var0, varargin )
      if ischar(var0)
        var = strcmp(obj.Gvarnames,var0);
        assert( sum(var) == 1, sprintf('variable %s not found', var0) );
      else
        assert( var0 <= numel(obj.Gvarnames), 'invalid variable number' );
        var = var0;
      end
      obj.G.imagesc( obj.GX(var), varargin{:} );
      title( obj.Gvarnames{var} );
    end
    
    function plotvars( obj )
      figI = 1;
      for i = 1:numel(obj.Gvarnames)
        if strcmp('intercept',obj.Gvarnames{i}); continue; end
        figure(figI); clf; obj.plotvar(i);
        figI = figI + 1;
      end
    end
    
    function save_variable_images( obj )
      for i = 1:numel(obj.Gvarnames)
        vname = obj.Gvarnames{i};
        if strcmp('intercept',vname); continue; end
        obj.G.save_image( obj.GX(i), [class(obj), '-variable'], 'name', vname );
      end
    end
    
  end
  
end
