function Z = project( obj, X0, Y0, Z0 )
% project data Z with coordinates X,Y onto the grid

assert(isrow(X0));
assert(isrow(Y0));
assert(isrow(Z0));

% inside grid?
assert( min(X0) >= obj.xgMin );
assert( max(X0) <= obj.xgMax );
assert( min(Y0) >= obj.ygMin );
assert( max(Y0) <= obj.ygMax );

sumZ0 = sum(Z0(:));

[I,J] = obj.coord2sub( X0, Y0 );
Z = zeros(size(obj));

% fix grid borders
sel = Z0 > 0 & ( rem(I,1)>0 | rem(J,1)>0 );
if nnz(sel) > 0
  for el = find(sel)
    Z( floor(I(el)), floor(J(el)) ) = Z( floor(I(el)), floor(J(el)) ) + Z0(el) / 4;
    Z( floor(I(el)), ceil(J(el))  ) = Z( floor(I(el)), ceil(J(el))  ) + Z0(el) / 4;
    Z( ceil(I(el)) , floor(J(el)) ) = Z( ceil(I(el)) , floor(J(el)) ) + Z0(el) / 4;
    Z( ceil(I(el)) , ceil(J(el))  ) = Z( ceil(I(el)) , ceil(J(el))  ) + Z0(el) / 4;
    Z0(el) = 0;
  end
end

% split data on grid cells
sel = Z0 > 0;
for el = find(sel)
  Z(I(el), J(el)) = Z(I(el), J(el)) + Z0(el);
end

% sanity check
assert_zero_with_tolerance( sum(Z(:)) - sumZ0, 1e-8 );
end
