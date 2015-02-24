classdef Scattered < Domain.Domain
  %SCATTEREDDOMAIN represents data given as (X,Y,Z) tuples.
  %  Uses ScatteredInterpolant.
  
  properties
    datafile
    no_smoothing_method
  end
  
  methods
    function obj = Scattered( datafile )
      assert( ischar(datafile) );
      if ~exist(datafile,'file'); error('file doest not exist: %s', datafile); end
      obj.datafile = datafile;
      obj.no_smoothing_method = 'nearest';
    end
    
    function data = gridded( obj, variable, grid )
      obj.gridded@Domain.Domain(variable, grid);
      
      [xx,yy]=ndgrid( grid.xc, grid.yc );
      
      % data that are loaded via load_dataset
      if ~isa(variable.transform,'function_handle'); variable.transform = @sqrt; end
      switch variable.name
        case 'pop'
          dataset_name = variable.name;
        case 'psy'
          dataset_name = 'hosp';
      end
      if variable.smoothing
        switch variable.name
          case 'pop'
            data = load_dataset( dataset_name, grid, 'bw', variable.smoothing, 'kernel_function', variable.smoothingK, 'transform', variable.transform, 'offset', variable.offset );
            datafield_name = sprintf('dens%stot', func2str(variable.transform));
          case 'psy'
            data = load_dataset( dataset_name, grid, 'bw', variable.smoothing, 'kernel_function', variable.smoothingK, 'transform', variable.transform, 'offset', variable.offset );
            datafield_name = sprintf('dens%spsybeds', func2str(variable.transform));
          otherwise
            data = obj.readdata;
            datafield_name = variable.name;
        end
        data = data.(datafield_name);
      else
        % no smoothing
        switch variable.name
          case 'pop'
            data = importdata(obj.datafile);
            datafield_name = 'tot';
          case 'psy'
            data = importdata(obj.datafile);
            datafield_name = 'psybeds';
          otherwise
            data = obj.readdata;
            datafield_name = variable.name;
        end
        z = double( data.(datafield_name) );
        % apply transform before interpolation
        if isa(variable.transform,'function_handle')
          z = variable.transform(z + variable.offset);
        end
        SI = scatteredInterpolant( data.x', data.y', z', obj.no_smoothing_method, 'none' );
        data = SI(xx,yy);
      end
      
    end
    
    function Vq = value( obj, variable, Xq, Yq )
      obj.value@Domain.Domain( variable, Xq, Yq );
      % data that are loaded via load_dataset
      if ~isa(variable.transform,'function_handle'); variable.transform = @sqrt; end
      switch variable.name
        case 'pop'
          Vq = importdata(obj.datafile);
          datafield_name = 'tot';
        case 'psy'
          Vq = importdata(obj.datafile);
          datafield_name = 'psybeds';
        otherwise
          Vq = obj.readdata;
          datafield_name = variable.name;
      end
      z = double( Vq.(datafield_name) );
      if isa(variable.transform,'function_handle')
        % apply transform before smoothing
        z = variable.transform(z + variable.offset);
      end
      if variable.smoothing
        Vq = nadaraya_watson_interpolation( Vq.x', Vq.y', z', Xq, Yq, variable.smoothing );
      else
        % no smoothing
        if variable.extrapolate; extr = 'nearest'; else extr = 'none'; end
        SI = scatteredInterpolant( Vq.x', Vq.y', z', 'nearest', extr );
        Vq = SI(Xq,Yq);
      end
      
    end
    
    function bb = bounding_box( obj, variable )
      data = importdata(obj.datafile);
      bb = [min(data.x) min(data.y); max(data.x) max(data.y)];
    end
    
    function data = readdata( obj )
      %% read datafile
      NANVALUE = int32(-2147483648);
      data = load(obj.datafile);
      data = parse_data( data );
      
      function data = parse_data( data )
        for f = fieldnames(data)'
          f = f{:};
          coldat = data.(f);
          colclass = class(coldat);
          switch colclass
            case 'int32'
              nans = coldat==NANVALUE;
              data.(f) = double( data.(f) );
              data.(f)(nans) = nan;
            case 'double'
            case 'struct'
              data.(f) = parse_data( coldat );
            otherwise
              error('Unsupported column format "%s" for column "%s"', colclass, f);
          end
        end
      end
    end
  end
  
end
