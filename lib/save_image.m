function save_image( dens, name, varargin )
%save_image Summary of this function goes here
%   Detailed explanation goes here

% check input structure
if isstruct(dens)
  x = dens.x;
  y = dens.y;
  dens = dens.z;
  assert( numel(dens) == numel(x)*numel(y) );
end


if isodd(nargin)
  error('Variable number of arguments must be even.');
end

Ncolors = 256;

%% build output filename
imwrite_name = name;
for k=1:2:length(varargin)-1
  imwrite_name = [imwrite_name, '-', varargin{k}, '=', num2str( varargin{k+1},'%0.3g' )];
  if strcmp( varargin{k}, 'eps' ) || strcmp( varargin{k}, 'grid_width' )
    Eps = varargin{k+1};
  end
end

%% calculate image
value_range_loc = strcmp('value-range',varargin);
if any(value_range_loc)
  value_range = varargin{find(value_range_loc)+1};
  assert( all(size(value_range)==[1 2]), 'Need a min. and a max. value' );
  interval_len = range(value_range);
  interval_min = min(value_range);
  interval_max = max(value_range);
  if max(dens(:)) > interval_max || min(dens(:)) < interval_min
    warning('Values outside predefined range.');
  end
  dens = round( flip( dens )/ interval_len *(Ncolors-1) );
  dens( dens < 0 ) = 0;
  dens( dens > Ncolors ) = Ncolors;
else
  interval_len = len(dens(dens(:)<Inf & dens(:)>-Inf));
  if isempty(interval_len) || interval_len==0
    warning('density %s is invalid',name);
  else
    interval_min = min(interval_len);
    dens = round( flip( dens )/ interval_len *(Ncolors-1) );
    dens = dens - min(dens(dens(:)>-Inf));
  end
end

%% negative Werte können Vorkommen
if any( isnan(dens(:)) )
  disp(['Warning: have NaNs; mapping to color 0 - ', imwrite_name]);
  dens(isnan(dens)) = 0;
end

if any(isinf(dens(:)))
  disp(['Warning: have Inf; mapping to min/max. color - ', imwrite_name]);
  dens(dens==Inf) = Ncolors - 1;
  dens(dens==-Inf) = 0;
end

assert( all(dens(:) >= 0 ) );
assert( all(dens(:) <= Ncolors-1 ) );
newmap = ( jet(Ncolors) );

%% fix filename for saving
% replace spaces
imwrite_name = regexprep(imwrite_name, '\s+', '_');
imwrite_args = [imwrite_name, '.png'];

%% save image
imwrite(dens, newmap, imwrite_args);

%% save metadata to file
if exist('Eps','var') && exist('x','var')
  fid = fopen([imwrite_name, '.pgw'], 'wt' );
  fprintf( fid, '%f\n%f\n%f\n%f\n%f\n%f\n', Eps*1000, 0, 0, -Eps*1000, min(x)*1000, max(y)*1000 );
  fclose(fid);
end
end
