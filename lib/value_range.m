function range = value_range( Cell )
assert( iscell( Cell ), 'Need a cell array.' );
value_range_min = min( cellfun( @(x) min(x(:)), Cell(:) ) );
value_range_max = max( cellfun( @(x) max(x(:)), Cell(:) ) );
range = [value_range_min, value_range_max];
end
