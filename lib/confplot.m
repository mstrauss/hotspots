function [ P ] = confplot( X, Y, YE, with_smoothing )
assert( all( size(X) == size(Y) ) );
assert( all( size(X) == size(YE) ) );

Y0 = Y;
YE = YE*1.96;
Yupper = Y+YE;
Ylower = Y-YE;
if ~exist('with_smoothing','var'); with_smoothing = true; end
if with_smoothing
  Y = smooth(Y)';
  Yupper = smooth(Yupper)';
  Ylower = smooth(Ylower)';
end

fill_between_lines = @(X,Y1,Y2,C) fill( [X fliplr(X)],  [Y1 fliplr(Y2)], C );
fill_between_lines( X, Yupper, Ylower, [0.7 0.7 0.7] ); hold on
P = plot( X, Y0, 'x' );
legend(P, 'data');
if with_smoothing
  p=plot(X,Y,'b:');
  legend(p, 'smoothed data')
  P = [P p];
end
legend(P);
hold off
end
