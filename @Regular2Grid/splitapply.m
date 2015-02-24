function varargout = splitapply( grid, N, fun, maxProd )
%% splitapply splits the grid into blocks and applies fun to it
% N: number of (unsplittable) input records
% it returns an array of grid size with the results of fun
% fun must take a bounding_box_ind as single parameter

if ~exist('maxProd','var'); maxProd = 1e8; end

Nq = grid.Nx * grid.Ny;

if N*Nq > maxProd
  iterations = ceil( N*Nq / maxProd );
  if isodd(iterations); iterations = iterations + 1; end
else
  iterations = 1;
end

if iterations > 1
  xFactors = factor(grid.Nx);
  Ncols = xFactors(1);
  Nrows = iterations/Ncols;
  yFactors = factor(grid.Ny);
  Nrows = yFactors( find( yFactors-Nrows > 0, 1 ) );
  fprintf('Regular2Grid#splitapply: Splitting [%d x %d] array into %dx%d pieces.\n', N, Nq, Nrows,Ncols);
else
  Ncols = 1;
  Nrows = 1;
end

for i = 1:nargout
  varargout{i} = [];
end

for column = 1:Ncols
  for i = 1:nargout
    Vqcol{i} = [];
  end
  for row = 1:Nrows
    bbi = [...
      (column-1)*grid.Nx/Ncols+1, ...
      (row-1)*grid.Ny/Nrows+1; ...
      column*grid.Nx/Ncols, ...
      row*grid.Ny/Nrows];
    res = fun(bbi);
    for i = 1:nargout;  Vqcol{i} = [Vqcol{i} squeeze(res(i,:,:))];  end
    fprintf( '%2dx%2d √\n', row, column );
    %     for i = 1:nargout
    %       figure(i)
    %       imagesc(Vqcol{i}'); axis xy; colorbar; drawnow;
    %     end
  end
  for i = 1:nargout;  varargout{i} = [varargout{i}; Vqcol{i}];  end
end

for i = 1:nargout
  assert( all( size(varargout{i}) == size(grid) ) );
end
