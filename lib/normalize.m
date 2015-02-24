function [ X, N ] = normalize( X )
%NORMALIZE a matrix of row vectors
N = sqrt(sum(X.^2,2));
N(N==0)=1;
X = bsxfun( @rdivide, X, N );
end
