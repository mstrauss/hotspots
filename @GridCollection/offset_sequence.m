function [seq, bfx, bfy] = offset_sequence( Ngrids )
%% static class method
% seq: sequence of grid offsets
% bfx, bfy: blowup factors

switch Ngrids
  case 5
    seq = [0 0; 2 0; 1 1; 0 2; 2 2];
    bfx = 3;
%   case 9
%     seq = [0 0; 0 2; 3 1; 0 2; 2 2; 4 2; 1 3; 3 3; 2 4];
%     bfx = 5;
  otherwise
    bfx = sqrt(Ngrids);
    if mod( bfx, 1 ) == 0
      seq = combinations(0:bfx-1, 0:bfx-1);
    else
      error('Unsupported number of grids: %d', Ngrids);
    end
end

bfy = bfx;
