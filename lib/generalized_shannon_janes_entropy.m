function [S,Ssum] = generalized_shannon_janes_entropy( A, M )
% % fix data
% A = A ./ sum(A(:));
% M = M ./ sum(M(:));
% assert( all( A(:) >= 0 ) );
% assert( all( M(:) >= 0 ) );

S = A - M - A .* log( A ./ M );
S( isnan(S) ) = 0;
S( isinf(S) ) = 0;
%infS = find(isinf(S(:)));
%if(infS) warning('Infinity at %d\n', infS); end
Ssum = sum(S(:));
end
