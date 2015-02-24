function marked_network = wiggled_marks_random_field( marked_network, wiggle_amount )
assert( isa(marked_network, 'MarkedNetwork') );
wiggle_ed = marked_network.net.rmove( marked_network.edge_network, wiggle_amount );
marked_network = MarkedNetwork( marked_network.net, wiggle_ed );
end
