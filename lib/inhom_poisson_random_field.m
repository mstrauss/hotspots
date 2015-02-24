function marked_network = inhom_poisson_random_field( network, number_of_points, grid, intensity )
assert( isa(network, 'Network') );
assert( isa(grid, 'Regular2Grid') );
[~,~,marked_network] = network.rpoisinhom( number_of_points, intensity, grid );
end
