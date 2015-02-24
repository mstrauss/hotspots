function [ norms ] = vnorm( X )
%vnorm Calculate norm of matrix of row vectors.

norms = sqrt(sum(X.^2,2));

end
