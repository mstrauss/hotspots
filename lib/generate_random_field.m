function [gHistory, dHistory, nnHistory, fields] = generate_random_field( number_of_fields, generating_function, ctr )
% gHistory: g(r), pair distribution function
% dHistory: d(r), nearest neighbor distr. fct.
% nnHistory: nn(k), mean distance to k-th neighbor

%assert( number_of_fields >= 39 );

% iteration count
N = 0;
m = numel(ctr);
gHistory = zeros(number_of_fields,m);
dHistory = zeros(number_of_fields,m);
nnHistory = zeros(number_of_fields,m);

if nargout == 4
  fields = cell(1,number_of_fields);
end

% generate fields & collect statistics
for i = 1:number_of_fields
  field = generating_function(i);
  if nargout == 4
    fields{i} = field;
  end
  % g(r)
  g = pcf(field, ctr);
  gHistory(i,:) = g;
  % d(r)
  [n1,~,s] = field.nearest_neighbor_distances;
  d = hist( n1, ctr );
  dHistory(i,:) = d;
  % nn(k)
  s(isinf(s))=0;
  nn = sum(s(:,1:m))/(size(s,1)-1);
  nnHistory(i,:) = nn;
  
  N = N + 1;
  % indicate progress
  if ~mod(i,round(number_of_fields/20)); fprintf('%.0f%% ', i/number_of_fields*100); end
end
fprintf('\n');

end
