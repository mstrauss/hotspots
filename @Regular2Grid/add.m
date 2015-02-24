function [ resulting_data, resulting_grid ] = add( obj, other_grid, data1, data2 )

assert( numel(data1) == obj.Nx*obj.Ny );
assert( numel(data2) == other_grid.Nx*other_grid.Ny );

[ resulting_grid, bf1x, bf1y, bf2x, bf2y] = obj.combine( other_grid );

m1 = kron( data1, repmat( 1, bf1x, bf1y ) );
m2 = kron( data2, repmat( 1, bf2x, bf2y ) );

m1 = padarray( m1, positive([resulting_grid.Nx resulting_grid.Ny]-size(m1))/2 );
m2 = padarray( m2, positive([resulting_grid.Nx resulting_grid.Ny]-size(m2))/2 );

resulting_data = m1 + m2;

  function x = positive( x )
    x(x<0) = 0;
  end
end

