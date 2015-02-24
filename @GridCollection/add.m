function [ resulting_data, dim3_sum ] = add( obj, data )
% add the data for the grids (which have been given at object construction)

warning('Regular2Grid#add is obsolete. Use "quantile" instead.');

prepared_data = obj.prepare( data );
resulting_data = sparse(obj.G.Nx, obj.G.Ny);
[m,n] = size(prepared_data);
dim3_sum = zeros(obj.G.Nx, obj.G.Ny, m);
for i = 1:m  % iterate over records
  for j = 1:n  % iterate over grids
    dim3_sum(:,:,i) = dim3_sum(:,:,i) + prepared_data{i,j};
    resulting_data = resulting_data + prepared_data{i,j};
  end
end

assert( all( size(resulting_data) == size(obj.G) ) );

end
