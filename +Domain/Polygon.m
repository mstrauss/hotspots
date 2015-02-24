classdef Polygon < Domain.Domain
  
  properties %(Access=private)
    shapefile
    edge_value_cache
    edge_value_cache_ref
  end
  
  methods
    function obj = Polygon( shapefile )
      if ~exist(shapefile,'file')
        error('shapefile "%s" not found',shapefile);
      end
      obj.shapefile = shapefile;
    end
    
    function data = gridded( obj, variable, grid )
      data = obj.gridded@Domain.Domain(variable, grid);
      
      [xx,yy]=ndgrid( grid.xc, grid.yc );
              
      [poly, map] = grid.find_shape_polygons(obj.shapefile);
      if variable.categorical
        values = categorical(data);
        values(poly>0)=categorical({map(poly(poly>0)).(variable.name)});
        cats = categories(values);
        cats = cats(2:end);
        data = cell(numel(cats),1);
        for i = 1:numel(cats)
          data{i} = struct();
          data{i}.data = zeros(size(grid));
          level = cats{i};
          data{i}.data = double(values==level);
          data{i}.level = level;
        end
        if isa(variable.transform,'function_handle')
          warning('Variable %s has a transform but is marked categorical. Ignoring transform.');
        end
      elseif isa(variable.transform,'function_handle')
        data(poly>0) = [map(poly(poly>0)).(variable.name)];
        data = variable.transform(data + variable.offset);
      end 
           
    end
    
    function [Vq, poly, map] = value( obj, variable, Xq, Yq, poly, map )
      obj.value@Domain.Domain( variable, Xq, Yq );
      if ~exist('poly','var') || ~exist('map','var')
        [poly, map] = obj.find_shape_polygons(Xq, Yq);
      end
      if variable.categorical
        values = categorical( {map.(variable.name)} );
        Vq = categorical(zeros(numel(Xq),1));
        Vq(poly>0) = values(poly(poly>0));
      else
        Vq = zeros(numel(Xq),1);
        Vq(poly>0) = [map(poly(poly>0)).(variable.name)];
        if isa(variable.transform,'function_handle')
          Vq = variable.transform(Vq + variable.offset);
        end
      end
      if variable.extrapolate
        sel = poly==0;
        if nnz(sel) > 0
          Vq(sel) = obj.extrapolate( map, Xq(sel), Yq(sel), variable );
          if isa(variable.transform,'function_handle')
            Vq(sel) = variable.transform(Vq(sel) + variable.offset);
          end
        end
      end
    end
    
    function Vq = edge_value( obj, variable, edge_coordinates )
      obj.edge_value@Domain.Domain( variable, edge_coordinates )
      if edge_coordinates == obj.edge_value_cache_ref
        % load from cache
        poly = obj.edge_value_cache.poly;
        map = obj.edge_value_cache.map;
        Xq = obj.edge_value_cache.Xq;
        Yq = obj.edge_value_cache.Yq;
        Vq = obj.value( variable, Xq, Yq, poly, map );
      else
        [Xq,Yq] = edge_coordinates.XY;
        [Vq, poly, map] = obj.value( variable, Xq, Yq );
        % update cache
        obj.edge_value_cache_ref = edge_coordinates;
        obj.edge_value_cache = struct('poly', poly, 'map', map, ...
          'Xq', Xq, 'Yq', Yq );
      end
    end
    
    function bb = bounding_box( obj, variable )
      S = obj.readdata;
      bb = [min([S.X]) min([S.Y]); max([S.X]) max([S.Y])];
    end
    
    function data = readdata( obj )
      data = shaperead(obj.shapefile);
    end
    
  end
  
  methods (Access=private)
    
    function [poly, map] = find_shape_polygons( obj, X, Y )
      map = obj.readdata;
      poly = map_assign_polygons( map, X, Y );
    end
    
    function extrapolated_data = extrapolate( ~, map, Xq, Yq, variable )
      % extrapolate by finding the closest polygon and taking its data
      mindist = Inf(1, numel(Xq));
      closestpoly = zeros(1, numel(Xq));
      for i = 1:numel(map)  % loop over polygons
        d = pdist2([map(i).X' map(i).Y'], [Xq Yq], 'euclidean', 'Smallest', 1);
        sel = d < mindist;
        mindist(sel) = d(sel);  closestpoly(sel) = i;
      end
      assert( all(closestpoly > 0), 'no closest polygon found!' );
      extrapolated_data = [map(closestpoly).(variable.name)]';
      assert( all(size(extrapolated_data) == size(Xq)) );
    end
  end
end
