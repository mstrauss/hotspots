function quan = quantile( obj, data, quantiles, bbi )
%% calculate the subpixel-wise lower q-quantile
% data: cell array of cell arrays containing the subgrids
% ind: linear index of grid cells to calculate

% cell array of prepared subpixel data
if ~exist('bbi','var'); bbi = obj.bounding_box_ind; end
prepared_data = obj.prepare( data, bbi );

[m,n] = size(prepared_data);
assert( n == obj.N, 'should equal number of grids' );

Nx = bbi(2)-bbi(1)+1;  Ny = bbi(4)-bbi(3)+1;  N = m*n;


M=reshape(full(cell2mat(prepared_data(:)')), Nx, Ny, N);
S = sort(M,3);
for i = 1:numel(quantiles);
  q = quantiles(i);
  if q==0
    quan(i,:,:) = mean(M,3);
  else
    qi = N*q;
    if mod(qi,1)
      qi = [floor(qi) ceil(qi)];
      if any(qi>0)
        warning('quantile cut off');
        qi = qi(qi>0);
      end
      quan(i,:,:) = mean( S(:,:,qi), 3 );
    else
      quan(i,:,:) = S(:,:,qi);
    end
  end
end
