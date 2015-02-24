function save_figure( name, varargin )
%save_figure

if ~isodd(nargin)
  error('Variable number of arguments must be odd.');
end

imwrite_name = name;
for k=1:2:length(varargin)-1
  imwrite_name = [imwrite_name, '-', varargin{k}, '=', num2str( varargin{k+1},'%0.2g' )];
end
% replace spaces
imwrite_name = regexprep(imwrite_name, '\s+', '_');

imwrite_args = [imwrite_name, '.png'];

saveas(gcf, imwrite_args);
end
