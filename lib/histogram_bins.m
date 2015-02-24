function ctr = histogram_bins( binwidth, maximum )
ctr = binwidth/2:binwidth:maximum;
assert( ctr(1) == (ctr(2)-ctr(1))/2 );
end
