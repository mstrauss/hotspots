function [Poly, S] = find_shape_polygons( grid, shapefile )
S = shaperead(shapefile);

% load from cache
global cachedir
hash = sum(shapefile) + grid.Nx + grid.Ny;
cachefile = fullfile( cachedir, sprintf('find_shape_polygons_cache_%d.mat',hash));
if exist(cachefile,'file')
  load( cachefile )
end

if ~exist('Poly','var')
  fprintf('Calculating polygons (nothing cached, hash=%d).\n', hash);
  recalc_poly
elseif numel(Poly) ~= grid.Nx * grid.Ny
  fprintf('Recalculating polygons (cache [hash=%d] mismatch).\n', hash);
  recalc_poly
end

  function recalc_poly
    [X,Y] = ndgrid( grid.xc, grid.yc );
    Poly = map_assign_polygons( S, X, Y );
    % update cache
    save(cachefile,'Poly');
  end

end
