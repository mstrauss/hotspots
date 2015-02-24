function [ n, ctr, Y ] = pcf( data, ctr, nbins )
% PCF calculates the pair distribution from a distance matrix D.
%   Return values:
%     n: counts
%     ctr: bin centers
%     Y: raw distance vector from which the histogram is made
%        (can be used in ksdensity to get a smoother pcf)

if ~exist('nbins','var')
  nbins = 500;  % 200 bins should be enough for smooth curve
end

if exist('ctr','var')
  assert_zero_with_tolerance( ctr(1) - (ctr(2)-ctr(1))/2, 1e-8 );
end

if isa(data,'MarkedNetwork')
  D = data.MD;
  assert( data.is_simple );
  %F = data.freq;
  F = ones( size(D,1), 1 );
else
  D = data;
  F = ones( size(data,1), 1 );
end

N = size(D,1);
assert( N == size(D,2) );  % D must be a square matrix
U = triu(D,1);
% fix potential multiplicity of points
I = repelem( (1:N)', F );
U = U(I,I);

Y = U(:);
% remove Inf and Zeros
Y = Y( Y>0 & Y<Inf );
if ~exist('ctr','var')
  binwidth = max(Y(:)) / nbins
  ctr = binwidth/2:binwidth:max(Y(:))
end
[n, ctr] = hist( Y, ctr );
%semilogy(ctr, n./ctr, 'x');
end
