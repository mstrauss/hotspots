function [ ] = save_vars_text( filename, varargin )

if isodd(length(varargin))
  error('Variable number of arguments must be even.');
end

if ~isempty(varargin)
  fid = fopen(filename, 'wt' );
  for k=1:2:length(varargin)-1
    v = varargin{k+1};
    if iscell(v)
      v = cell2table(v);
    end
    if istable(v)
      writetable(v, sprintf('%s-%s.txt', filename, varargin{k}));
    elseif numel(v) > 1
      fprintf( fid, '%s = %s\n', varargin{k}, mat2str(v) );
    elseif isa(v,'function_handle')
      fprintf( fid, '%s = %s\n', varargin{k}, func2str(v) );
    else
      fprintf( fid, '%s = %.6g\n', varargin{k}, v );
    end
  end
  fclose(fid);
end


end

