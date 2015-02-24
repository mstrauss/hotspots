function [ new_outdir, dirname ] = mk_outdir( varargin )
project_configuration
if isodd(nargin)
  error('Variable number of arguments must be even.');
end

dirname = datestr(now,30);
if ~isempty(varargin)
  subname = '';
  for k=1:2:length(varargin)-1
    vardat = varargin{k+1};
    if isa(vardat,'function_handle')
      strvardat = func2str(vardat);
    else
      strvardat = num2str( vardat,'%0.3g' );
    end
    subname = [subname, '-', varargin{k}, '=', strvardat];
  end
  % replace spaces
  subname = regexprep(subname, '\s+', '_');
  dirname = sprintf('%s%s', dirname, subname);
end

if strcmp(outdir,fullfile(pwd,'results/'))
  new_outdir = fullfile( outdir,dirname );
  mkdir( new_outdir );
  fprintf('setting outdir = %s\n', new_outdir);
end

end
