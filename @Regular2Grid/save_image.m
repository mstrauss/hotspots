function save_image( obj, data, name, varargin )
if size(data')==size(obj); data = data'; end
assert( all(size(data)==size(obj)) );
save_image( struct('x', obj.xc, 'y', obj.yc, 'z', data' ), ...
  name, 'grid_width', obj.ex, varargin{:} );

end
