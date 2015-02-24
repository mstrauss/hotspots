classdef Variable < handle
  %VARIABLE describes a network variable.
  
  properties
    name          % name of the variable
    domain        % the domain of the definition of the variable
    % (subclass of Domain)
    categorical   % true if categorical
    transform     % transformation function_handle (1 parameter) or false
    smoothing     % smoothing bandwidth or false
    smoothingK    % smoothing kernel function name
    offset        % an offset, being applied before transformation
    extrapolate   % true if extrapolation should be done
    fixNan        % true if NaN shall be fixed
  end
  
  methods (Static)
    function transforms = transforms()
      % valid transforms (or any function_handle with one parameter)
      transforms = {'log', 'sqrt'};
    end
  end
  
  methods
    
    function obj = Variable( name, varargin )
      
      obj.name = name;
      
      % parse options
      argi = 1;
      while argi <= numel(varargin)
        switch varargin{argi}
          case 'domain'
            argi = argi + 1; obj.domain = varargin{argi};
          case 'categorical'
            argi = argi + 1; obj.categorical = varargin{argi};
          case 'transform'
            argi = argi + 1; obj.transform = varargin{argi};
          case 'smoothing'
            argi = argi + 1; obj.smoothing = varargin{argi};
          case 'smoothingK'
            argi = argi + 1; obj.smoothingK = varargin{argi};
          case 'offset'
            argi = argi + 1; obj.offset = varargin{argi};
          case 'extrapolate'
            argi = argi + 1; obj.extrapolate = varargin{argi};
          case 'fixNan'
            argi = argi + 1; obj.fixNan = varargin{argi};
          otherwise
            error('Unkown option %s',varargin{argi});
        end
        argi = argi + 1;
      end
      
      % find a domain
      if isempty(obj.domain)
        switch obj.name
          case 'intercept'
            obj.domain = Domain.Constant;
          otherwise
            error('require domain for variable %s', obj.name);
        end
      end
      
      % default options
      if isempty(obj.categorical); obj.categorical = false; end
      if isempty(obj.transform); obj.transform = false;
      else
        assert( isa(obj.transform, 'function_handle') || ismember(obj.transform, Variable.transforms), 'invalid transform' );
        if ~isa(obj.transform, 'function_handle')
          obj.transform = str2func(obj.transform);
        end
      end
      if isempty(obj.smoothing); obj.smoothing = false; end
      if isempty(obj.smoothingK); obj.smoothingK = 'gauss'; end
      if isempty(obj.offset); obj.offset = 0; end
      if isempty(obj.extrapolate); obj.extrapolate = false; end
      if isempty(obj.fixNan); obj.fixNan = false; end
    end
    
    function data = gridded( obj, grid )
      %% delegate to domain
      data = obj.domain.gridded( obj, grid );
      data = obj.fixnan( data );
    end
    
    function data = edge_value( obj, edge_coordinates )
      %% find the value of the variable at given edge coordinates
      data = obj.domain.edge_value( obj, edge_coordinates );
      data = obj.fixnan( data );
      if ~iscolumn(data); data = data'; end
    end
    
    function Vq = value( obj, Xq, Yq )
      %% return the value of the variable at (Xq,Yq)
      Vq = obj.domain.value( obj, Xq, Yq );
      Vq = obj.fixnan( Vq );
    end
    
    function bb = bounding_box( obj )
      %% returns the bounding box of this variable
      bb = obj.domain.bounding_box( obj );
    end
  end
  
  methods (Access=private)
    function data = value_on_mark( obj, mark_index )
      %% retrieve the value when it is defined for a mark
      G = obj.domain.graph;
    end
    
    function data = value_on_edge( obj, edge_index )
      %% retrieve the value when it is defined for an edge
      net = obj.domain.net;
      net.edgedata.(obj.name)
    end
    
    function data = fixnan( obj, data )
      if obj.fixNan
        replacement = 0;
        if isa(obj.transform,'function_handle')
          replacement = obj.transform(replacement + obj.offset);
        end
        nans = isnan(data);
        if nnz(nans)>0
          fprintf('Variable "%s": Fixing NaNs as requested by %f.\n', obj.name, replacement);
          data(nans) = replacement;
        end
      end
    end
  end
  
end

