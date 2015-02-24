function [ data, G ] = load_dataset( type, grid, varargin )
% load_dataset load a precalculated dataset from disk or calucates the
% required data.

project_configuration

if nargin < 2
  data = load_default_data(type);
elseif nargin == 2
  % only load default data, no densities, and a grid
  if isa(grid,'Regular2Grid')
    assert( grid.ex == grid.ey, 'require quadratic grid' );
    grid = grid.ex;
  end
  [data, G] = load_default_data(type, grid);
else
  % if we need more, we have have to specify a grid
  assert( isa(grid,'Regular2Grid'), 'require grid as second parameter' );
  i = 1;
  while i <= length(varargin)
    switch(varargin{i})
      case 'kernel_function'
        i=i+1; kernel_function = varargin{i};
      case 'force_calc'
        force_calc = true;
      case 'bw'
        i=i+1; bw = varargin{i};
      case 'transform'
        i=i+1; transform = varargin{i};
      case 'offset'
        i=i+1; offset = varargin{i};
      otherwise
        error('Unknown parameter %s', varargin{i});
    end
    i=i+1;
  end
  
  if ~exist('kernel_function', 'var'); kernel_function = 'gauss'; end
  if ~exist('force_calc', 'var'); force_calc = false; end
  if ~exist('transform', 'var'); transform = ''; end
  if ~exist('offset', 'var'); offset = 0; end
  
  data_file = data_filename( type, grid.ex, bw, kernel_function, transform );
  
  if ~force_calc && exist(data_file, 'file')
    disp(['Loading file ',data_file]);
    data=importdata( data_file );
    if ~isfield(data, 'kernel'); data.kernel = 'gauss'; end
    if ~isfield(data, 'offset'); data.offset = 0; end
    try
      assert( strcmp( data.kernel, kernel_function ), 'kernel mismatch' );
      assert( ~isfield(data,'transform') || strcmp(data.transform, func2str(transform)), 'invalid transform' );
      assert( all( [data.Nx, data.Ny] == size(grid) ), 'incompatible grid' );
      assert( all( data.gx == grid.xc ), 'incompatible grid gx ≠ xc' );
      assert( all( data.gy == grid.yc ), 'incompatible grid gy ≠ yc' );
      assert( data.offset == offset, 'incompatible offset' );
    catch err
      warning('problem "%s" loading from file => recalculating', err.message);
      v = [ varargin 'force_calc' ];
      data = load_dataset( type, grid, v{:} );
      return
    end
      
  else
    if ~force_calc; disp(['Warning: file ', data_file, ' does not exist.']); end
    data = calc_data( type, grid, bw, kernel_function, transform, offset );
  end
end
end

function [data, G] = load_default_data( type, grid_width )
if exist('grid_width','var')
  assert( isnumeric( grid_width ), 'require numeric grid_width' );
end
global pop_epsg3035_default_file cases_epsg3035_default_file hospitals_default_file
switch type
  case 'pop'
    data = importdata(pop_epsg3035_default_file);
    data.tot = double(data.tot);
  case 'cases'
    data = importdata(cases_epsg3035_default_file);
    data.coords = [data.x' data.y'];
    data.N = numel(data.x);
  case 'hosp'
    data = importdata(hospitals_default_file);
    % clean data
    data.psybeds(isnan(data.psybeds))=0;
    data.nonpsybeds(isnan(data.nonpsybeds))=0;
  otherwise
    error('Unkown type: %s', type);
end

if exist('grid_width','var')
  % build grid / FIXME: hardcoded stuff
  lx = round(len(data.x)*0.05);
  ly = round(len(data.y)*0.05);
  G = Regular2Grid( grid_width, grid_width, min(data.x)-lx, min(data.y)-ly, max(data.x)+lx, max(data.y)+ly );
  data = annotate( data, G );
end
end

function data = calc_data( type, grid, bw, kernel_function, transform, offset )
fprintf('Calculating %s density for eps=%s, bw=%s', type, num2str(grid.ex), num2str(bw));
data = load_default_data( type, grid.ex );
data = annotate( data, grid );
switch type
  case 'pop'
    fields = {'tot'};
  case 'hosp'
    fields = {'psybeds' 'nonpsybeds'};
  case 'cases'
    [data.C,data.I,data.J,data.CM] = grid.boxcount(data.x,data.y);
    data.z = ones(1,data.N);
    fields = {'z'};
  otherwise
    error('Unknown type: %s', type);
end

if ~exist('transform','var'); transform = ''; end
if ~isempty(transform)
  if ischar(transform)
    transstr = transform;
    transform = str2func(transform);
  else
    transstr = func2str(transform);
  end
  for i = 1:numel(fields)
    % TRANSFORMING DATA here:
    data.([transstr,fields{i}]) = transform(offset + data.(fields{i}));
    % DO NOT fix infinities
    %     if any(isinf(data.([transform,fields{i}])(:)))
    %       warning('fixing infintities in field %s', [transform,fields{i}]);
    %       data.([transform,fields{i}])(isinf(data.([transform,fields{i}]))) = 0;
    %     end
  end
end

for i = 1:numel(fields)
  field = fields{i};
  if strcmp( kernel_function, 'gauss_convolute' )
    [data.(['dens',transstr,field]), data.(['dens',transstr,field,'_kernel'])] = grid.density( bw*eye(2), grid.project( data.x, data.y, data.([transstr,field]) ) );
  else
    [data.(['dens',transstr,field])] = gaussian_smoothing( data.x, data.y, grid.ex, kernel_function, bw, data.([transstr,field]), grid.xc, grid.yc )';
  end
end

data.bw = bw;
data.kernel = kernel_function;
data.offset = offset;

save( data_filename(type, grid.ex, bw, kernel_function, transform), 'data' );
end

function [ data_filename ] = data_filename( type, grid_width, bw, kernel_function, transform )
global cachedir
switch type
  case 'pop'
    prefix = 'population';
  case 'hosp'
    prefix = 'hospitals';
  otherwise
    prefix = type;
end

switch kernel_function
  case 'exp'
    kernel = '-kernel=exp';
  case 'box'
    kernel = '-kernel=box';
  case 'gauss_convolute'
    kernel = '-kernel=conv';
  otherwise
    kernel = '';
end

if ~exist('transform','var'); transform = ''; end
if ~isempty(transform)
  transform = ['-transform=', func2str(transform)];
end

data_filename = fullfile(cachedir, [prefix, '-eps=', num2str(grid_width), '-', num2str(bw, '%0.1f'), kernel, transform, '.mat']);
end

function data = annotate( data, grid )
data.Nx = grid.Nx;
data.Ny = grid.Ny;
data.gx = grid.xc;
data.gy = grid.yc;
data.ex = grid.ex;
data.ey = grid.ey;
data.N = numel(data.x);
end
