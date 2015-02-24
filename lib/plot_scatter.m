function [ output_args ] = plot_scatter( M, labels, Mtrue, Mtruetext, varargin )
% create a scatter plot of given data matrix using given labels
% if Mtrue is given, the true positions are marked

% parse options
argi = 1;
while argi <= numel(varargin)
  switch varargin{argi}
    case 'group'
      argi = argi + 1; group = varargin{argi};
    case 'Y'
      argi = argi + 1; Y = varargin{argi};
    case 'Ylabels'
      argi = argi + 1; Ylabels = varargin{argi};
    case 'title'
      argi = argi + 1; title_text = varargin{argi};
    otherwise
      error('Unkown option %s',varargin{argi});
  end
  argi = argi + 1;
end

% default options
if ~exist('group','var'); group = false; end
if ~exist('Y','var')
  Y = false;
  Ylabels = labels;
end
if ~exist('title_text','var'); title_text = 'Scatter Plot of Simulation Results'; end
if group
  if ~Y; Y = M; end
  [H,AX,BigAx] = gplotmatrix( M, Y, group, [], [], [], [], [], labels, Ylabels );
else
  [H,AX,BigAx] = plotmatrix( M );

  for i = 1:numel(labels)
    ylabel(AX(i,1), labels(i));
    xlabel(AX(end,i), labels(i));
  end
end

set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto');

if exist('Mtrue','var')
  for i = 1:numel(Mtrue)
    for j = 1:numel(Mtrue)
      if i==j continue; end
      hold( AX(i,j) )
      plot( AX(i,j), Mtrue(j), Mtrue(i), 'Xr' )
    end
  end
  if exist('Mtruetext','var') && ~isempty(Mtruetext)
    title_text = [title_text, ' / red: ', Mtruetext];
  end
end

title(BigAx, title_text);

end

