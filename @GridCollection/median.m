function [ mean_of_median, median_per_record ] = median( obj, data )
% add the data for the grids (which have been given at object construction)

warning('Regular2Grid#median is obsolete. Use "quantile" instead.');

prepared_data = obj.prepare( data );  % cell array of prepared subpixel data

[m,n] = size(prepared_data);
assert( n == obj.N, 'should equal number of grids' );

mean_of_median = zeros(obj.G.Nx, obj.G.Ny);
median_per_record = cell(m,1);
Nx = obj.G.Nx;  Ny = obj.G.Ny;
parfor ri = 1:m  % iterate over records
  % find median for each (sub)pixel from all grids
  M = reshape( full( cell2mat(prepared_data(ri,:)) ), Nx, Ny, n );
  assert( all( [size(M,1) size(M,2) size(M,3)] == [Nx Ny n] ) );
  m = median(M,3);
  mean_of_median = mean_of_median + m;
  median_per_record{ri} = sparse(m);
end

assert( all( size(mean_of_median) == size(obj.G) ) );
