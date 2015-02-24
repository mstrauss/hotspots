function marked_network = poisson_random_field( network, number_of_points )
assert( isa(network, 'Network') );
[~,marked_network] = network.rpois( number_of_points );
end
