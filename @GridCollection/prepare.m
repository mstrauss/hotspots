function [prepared_data] = prepare( obj, data, bbi, convolute )
% fix input cell array data data (for the grids which have been given at
% object construction) to have same size
%   bbi: bounding box indices of subgrid to prepare
%       [minXi minYi; maxXi maxYi]
%   convolute: boolean, convolute input data with weight matrix?
% Input format for data:
%   cell array of cell arrays;  outer cells hold the records, inner cells
%   the grids

N = obj.N;        % number of grids
M = numel(data);  % number of records

% input validation
assert(iscell(data));
for i = 1:M
  assert(all(size(data{i}) == [1 N]));
end

Gbbi = obj.bounding_box_ind;
if ~exist('bbi','var'); bbi = Gbbi; end
assert( all( bbi(1,:) >= Gbbi(1,:) ), 'bbi outside bounds' );
assert( all( bbi(2,:) <= Gbbi(2,:) ), 'bbi outside bounds' );
xiMin = bbi(1,1);  yiMin = bbi(1,2);
xiMax = bbi(2,1);  yiMax = bbi(2,2);
xN = xiMax - xiMin + 1;
yN = yiMax - yiMin + 1;

if ~exist('convolute','var'); convolute = false; end

for i = 1:N
  for j = 1:M
    assert( all( size(data{j}{i}) == size(obj.Grids{i}) ) );
  end
end

% mangle data
prepared_data = cell(M,N);  % first the records, then the subgrids
[offsets, bfx, bfy] = GridCollection.offset_sequence(N);
for i = 1:N
  for j = 1:M
    prepared_data{j,i} = sparse(xN, yN);
    tmp = kron( data{j}{i}, ones( bfx, bfy ) );
    if convolute
      [~,K] = obj.weights(i);
      tmp = conv2( tmp, K, 'same' );
    end
    %fprintf('%d/%d: size(tmp) = [%d %d]\n', i, j, size(tmp,1), size(tmp,2) );

    dx0 = offsets(i,1); dy0 = offsets(i,2);
    %fprintf('%d/%d: offset = [%d %d]\n', i, j, dx0, dy0 );
  
    prepared_data{j,i} = tmp(xiMin+dx0:xiMax+dx0, yiMin+dy0:yiMax+dy0);
  
    assert( all( size(prepared_data{j,i}) == [xN yN] ) );
  end
end

end
