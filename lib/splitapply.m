function Vq = splitapply( N, Nq, fun, maxProd )
%%
% N: number of (unsplittable) input records
% Nq: number of (splittable) output records
% maxProd: splitting threshold for N*Nq
% fun: function to apply

if ~exist('maxProd','var'); maxProd = 1e8; end

if N*Nq > maxProd
  iterations = ceil( N*Nq / maxProd );
  fprintf('splitapply: Splitting [%d x %d] array into %d pieces.\n', N, Nq, iterations);
else
  iterations = 1;
end

blocksize = ceil(Nq/iterations);
togo = Nq;
Vq = zeros(Nq,1);
for i = 1:iterations
  % select query points
  sel = ( 1:min(togo,blocksize) ) + (i-1) * blocksize;
  nSel = numel(sel);

  Vq(sel) = fun( sel );
  
  % update loop vars
  togo = togo - nSel;
  
  % progress
  fprintf('%4d ',i);
  if ~mod(i,10); fprintf('\n'); end
end

assert( togo == 0 );
assert( length(Vq) == Nq );
