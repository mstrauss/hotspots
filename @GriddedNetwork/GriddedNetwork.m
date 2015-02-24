classdef GriddedNetwork < handle
  %GRIDDEDNETWORK represents a Network viewed on a Regular2Grid.
  
  properties (SetAccess=immutable)
    network     % the network
    grid_width  % the input grid width
    grid        % the resulting grid   (or the other way round)
  end
  
  methods
    function obj = GriddedNetwork( network, grid )
      assert( isa(network, 'Network') );
      obj.network = network;
      assert( isa(grid,'Regular2Grid') || isnumeric(grid) );
      if isnumeric(grid)
        bb = obj.network.bounding_box;
        obj.grid = Regular2Grid( grid, grid, bb(1), bb(3), bb(2), bb(4) );
        obj.grid_width = grid;
      else
        obj.grid = grid;
        assert( grid.ex == grid.ey );
        obj.grid_width = grid.ex;
      end
    end
    
    function plot_netlen( obj, varargin )
      % parse options
      argi = 1;
      while argi <= numel(varargin)
        switch varargin{argi}
          case 'plot_network'
            argi = argi + 1; plot_network = varargin{argi};
          otherwise
            error('Unkown option %s',varargin{argi});
        end
        argi = argi + 1;
      end
      % default options
      if ~exist('plot_network','var'); plot_network = true; end
      % do stuff
      obj.grid.imagesc( ...
        obj.grid.line_density( ...
        obj.network.edges ));
      if plot_network
        hold on; obj.network.plot;
      end
    end
        
  end
  
end
