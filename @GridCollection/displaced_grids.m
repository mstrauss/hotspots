function Grids = displaced_grids( base_grid, Ngrids )
%% static class method
% generate N displaced grids from base_grid

Grids = cell(1,Ngrids);

[seq, bfx, bfy] = GridCollection.offset_sequence( Ngrids );

G = base_grid;
for i = 1:Ngrids
  dx = seq(i,1);  dy = seq(i,2);
  Grids{i} = Regular2Grid(G.ex, G.ey, ...
    G.xcMin - dx*G.ex/bfx, G.ycMin - dy*G.ey/bfy, ...
    G.xcMax + dx*G.ex/bfx, G.ycMax + dy*G.ey/bfy );
end
end
