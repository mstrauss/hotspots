function Poly = map_assign_polygons( S, X, Y )
% assign polygons
Poly = zeros(size(X));   % to which polygon does X,Y belong

%%% Variant looping over Polygons
N = numel(S);  % number of polygons
for i = 1:N  % loop over polygons
  bb = S(i).BoundingBox;
  Xr = X(~Poly);  Yr = Y(~Poly);
  inbb = Xr >= bb(1,1) & Xr <= bb(2,1) & Yr >= bb(1,2) & Yr <= bb(2,2);
  inpoly = inpolygon(Xr(inbb),Yr(inbb),S(i).X,S(i).Y);
  %fprintf('Cases in poly %d/%d: %d (numel: %d)\n', i, N, nnz(inpoly), numel(inpoly));
  I = find(~Poly);
  I = I(inbb);
  I = I(inpoly);
  Poly( I ) = i;
end

%%% Parallel Variant
% xMin = cellfun( @(bb) bb(1,1), {S.BoundingBox} );
% xMax = cellfun( @(bb) bb(2,1), {S.BoundingBox} );
% yMin = cellfun( @(bb) bb(1,2), {S.BoundingBox} );
% yMax = cellfun( @(bb) bb(2,2), {S.BoundingBox} );
% N = numel(X);
% parfor j = 1:N  % loop over points
%   inbb = X(j) >= xMin & X(j) <= xMax & Y(j) >= yMin & Y(j) <= yMax;
%   nBB = nnz(inbb);
%   if nBB>0
%     bbi = find(inbb);
%     for p = 1:nBB
%       inpoly = inpolygon(X(j),Y(j),S(bbi(p)).X,S(bbi(p)).Y);
%       if inpoly
%         Poly(j) = bbi(p);
%         break;  % for loop
%       end
%     end
%   end
% end

