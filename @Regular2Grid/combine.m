function [ resulting_grid, bf1x, bf1y, bf2x, bf2y] = combine( obj, other_grid )

% das eine Gitter muss ganzzahlig in das andere passen
assert( mod( max(obj.ex, other_grid.ex), min(obj.ex, other_grid.ex) ) == 0 );
assert( mod( max(obj.ey, other_grid.ey), min(obj.ey, other_grid.ey) ) == 0 );

if( obj.ex == other_grid.ex ); ex = obj.ex / 2; else
  ex = min(obj.ex, other_grid.ex);
end
if( obj.ey == other_grid.ey ); ey = obj.ey / 2; else
  ey = min(obj.ey, other_grid.ey);
end

bf1x = obj.ex / ex;
bf1y = obj.ey / ey;
bf2x = other_grid.ex / ex;
bf2y = other_grid.ey / ey;

x0 = min( obj.xgMin, other_grid.xgMin ) + ex/2;
y0 = min( obj.ygMin, other_grid.ygMin ) + ey/2;

x1 = max( obj.xgMax, other_grid.xgMax ) - ex/2;
y1 = max( obj.ygMax, other_grid.ygMax ) - ey/2;

resulting_grid = Regular2Grid( ex, ey, x0, y0, x1, y1 );

end
