function Z = interpolate( this_grid, this_grid_index, other_grid, other_grid_data, interpolation_method )
%INTERPOLATE (other_grid_data given on other_grid) on this_grid_index
%points given by flat indexing on this_grid.

assert( all( size(this_grid) >= size(other_grid) ), ...
  'interpolation grid must have same or higher resolution');

if ~exist('interpolation_method','var'); interpolation_method = 'linear'; end


% grid cell centers where we wanna interpolate
coords = this_grid.ind2coord( this_grid_index );

% other grid sampling coordinates
[xx,yy] = ndgrid(other_grid.xc, other_grid.yc);
Z = interpn( xx, yy, other_grid_data, coords(:,1), coords(:,2), interpolation_method );
