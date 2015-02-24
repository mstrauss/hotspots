function [ resulting_data, prepared_data ] = weighted_mean( obj, data )
% add the data for the grids (which have been given at object construction)

warning('Regular2Grid#weighted_mean is obsolete. Use "quantile" instead.');

prepared_data = obj.prepare( data );

% reweight data
for i = 1:obj.N
  %prepared_data(:,:,i) = prepared_data(:,:,i) .* obj.weights(i);
  [~,K] = obj.weights(i);
  prepared_data(:,:,i) = conv2( prepared_data(:,:,i), K, 'same' );
end
resulting_data = mean( prepared_data, 3 );
assert( all( size(resulting_data) == size(obj.G) ) );

end
