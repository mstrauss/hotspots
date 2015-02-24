function [ resulting_grid, bfx, bfy] = combine( obj )
% combine Grids into one

% number of Grids to combine
N = obj.N;

% check pairwise fitting
pairs = combnk(1:N, 2);
for i = 1:size(pairs,1);
  g1 = obj.Grids{pairs(i,1)};  g2 = obj.Grids{pairs(i,2)};
  assert( mod( max(g1.ex, g2.ex), min(g1.ex, g2.ex) ) == 0 );
  assert( mod( max(g1.ey, g2.ey), min(g1.ey, g2.ey) ) == 0 );
end

ex = Inf;  eX = cellfun( @(g) g.ex, obj.Grids);
ey = Inf;  eY = cellfun( @(g) g.ey, obj.Grids);

[~,bfx,bfy] = GridCollection.offset_sequence(N);

if var(eX) < 1e-10; ex = g1.ex / bfx; else
  ex = min(eX);
  warning('only equal grid spacing supported; this might fail');
end
if var(eY) < 1e-10; ey = g1.ey / bfy; else
  ey = min(eY);
  warning('only equal grid spacing supported; this might fail');
end

% validate blow factors
for i = 1:N
  assert( obj.Grids{i}.ex/ex == bfx );
  assert( obj.Grids{i}.ey/ey == bfy );
end

% build combined grid
x0 = max( cellfun( @(g) g.xgMin, obj.Grids) ) + ex/2;
y0 = max( cellfun( @(g) g.ygMin, obj.Grids) ) + ey/2;

x1 = min( cellfun( @(g) g.xgMax, obj.Grids) ) - ex/2;
y1 = min( cellfun( @(g) g.ygMax, obj.Grids) ) - ey/2;

resulting_grid = Regular2Grid( ex, ey, x0, y0, x1, y1 );

end
