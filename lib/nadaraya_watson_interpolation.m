function [ Vq, N ] = nadaraya_watson_interpolation( X, Y, V, Xq, Yq, bandwidth )
%% NADARAYA_WATSON_INTERPOLATION performs Nadaraya-Watson interpolation
% Nadaraya Watson interpolation is when the input coordinates have
% uncertainty associated with them.  Given points (X,Y) with values V(X,Y)
% and given Gaussian kernel bandwidth.  Calculates the interpolated /
% smoothed surface at points (Xq,Yq) using direct Nadaraya-Watson
% interpolation.  For large bandwidth, the interpolated values approach the
% mean value of the input.  For low bandwidth, the interpolated values have
% the same value as the input points at the input coordinates.

N = length(X);
assert( N == length(Y) );
assert( N == length(V) );

% clean input
V(isnan(V)) = 0;

Nq = length(Xq);
assert( Nq == length(Yq) );

Vq = splitapply( N, Nq, @blockfun, 1e8 );


  function output = blockfun( sel )
    nSel = numel(sel);
    
    % calculate all pairwise distances between (X,Y) and (Xq,Yq)
    D = pdist2([X Y], [Xq(sel) Yq(sel)]);
    
    % fold with Gaussian kernel
    bw2 = 2*bandwidth^2;
    K = exp( -D.^2/bw2 );
    Knorm = repmat( sum(K,1), N, 1);  % Nadaraya-Watson factor
    Ksel = K>0 & Knorm > 0;
    K(Ksel) = K(Ksel) ./ Knorm(Ksel);
    output = sum( repmat(V,1,nSel) .* K, 1 );
  end

end
