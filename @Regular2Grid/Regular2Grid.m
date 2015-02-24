classdef Regular2Grid < handle
  % RegularGrid describes an evenly-spaced, 2-dim. Grid
  % All Coordinates are centered in grid boxes!
  
  properties (SetAccess = immutable)
    ex, ey  % x and y distance
    Nx, Ny  % number of grid boxes [aka pixels] for each direction
    xcMin, ycMin  % start coordinates (lower left point; grid box centered)
    % derived properties
    xc, yc  % grid box center coordinates
    xcMax, ycMax  % center coord. maxima
    xg, yg  % grid line coordinates
    xgMin, xgMax, ygMin, ygMax
  end
  
  methods
    function obj = Regular2Grid( ex, ey, x0, y0, x1, y1 )
      obj.ex = ex;
      obj.ey = ey;
      obj.xcMin = x0;
      obj.ycMin = y0;
      obj.Nx = round( (x1-x0)/ex ) + 1;
      obj.Ny = round( (y1-y0)/ey ) + 1;
      % calculated
      obj.xcMax = (obj.Nx-1)*ex + x0;
      obj.ycMax = (obj.Ny-1)*ey + y0;
      obj.xc = x0 + [0 cumsum( repmat(ex,1,obj.Nx-1) )];
      obj.yc = y0 + [0 cumsum( repmat(ey,1,obj.Ny-1) )];
      obj.xg = x0 - ex/2 + [0 cumsum( repmat(ex,1,obj.Nx) )];
      obj.yg = y0 - ex/2 + [0 cumsum( repmat(ey,1,obj.Ny) )];
      obj.xgMin = min(obj.xg); obj.xgMax = max(obj.xg);
      obj.ygMin = min(obj.yg); obj.ygMax = max(obj.yg);
    end
    
    function X = size( obj )
      X = [obj.Nx obj.Ny];
    end
    
    function X = intersectLines( obj, lines )
      % Calculate all intersection points of the grid with given lines.
      %   lines: rows of [a b c]
      %   [lines implicitly given by formula: a*x + b*y = c]
      X = [];
      for i = 1:size(lines,1)
        line = lines(i,:);
        % vertical grid
        A = [1 0; line(1:2)];
        C = [obj.ex; 0] * (1:obj.Nx-1);
        C(2,:) = line(3);
        X = [X A\C];
        % horizontal grid
        A = [0 1; line(1:2)];
        C = [obj.ey; 0] * (1:obj.Ny-1);
        C(2,:) = line(3);
        Y = A\C;
        X = [X Y];
      end
      % remove invalid data
      out_of_bounds = X(1,:) < obj.xgMin | X(2,:) < obj.ygMin | X(1,:) > obj.xgMax | X(2,:) > obj.ygMax;
      X(:,out_of_bounds) = [];
      X = X';
    end
    
    function X = intersectLineSegments( obj, lineSegments )
      % Calculate all intersection points of the grid with given lines.
      %   lineSegments: rows of [x y dx dy]
      % Output format: rows of [xCoord yCoord lineIndex]
      
      assert( size(lineSegments,2) == 4 );
      
      % convert line segment to (a,b,c) form
      x = lineSegments(:,1); y = lineSegments(:,2);
      dx = lineSegments(:,3); dy = lineSegments(:,4);
      c = 1;
      a = 1 ./ ( x - y .* dx ./ dy );
      b = 1 ./ ( y - x .* dy ./ dx );
      
      % number of line segments
      N = size( lineSegments, 1 );
      
      xL = min( x, x + dx );
      xU = max( x, x + dx );
      yL = min( y, y + dy );
      yU = max( y, y + dy );
      
      X = [];
      for i = 1:N
        % vertical grid
        xi = obj.xg >= xL(i) & obj.xg <= xU(i);
        %xi = 1:obj.Nx;
        A = [1 0; a(i) b(i)];
        C = [1; 0] * obj.xg(xi);
        C(2,:) = c;
        X = [X [A\C; repmat(i,1,size(C,2))]];
        % horizontal grid
        yi = obj.yg >= yL(i) & obj.yg <= yU(i);
        %yi = 1:obj.Ny;
        A = [0 1; a(i) b(i)];
        C = [1; 0] * obj.yg(yi);
        C(2,:) = c;
        X = [X [A\C; repmat(i,1,size(C,2))]];
      end
      % remove invalid data
      out_of_bounds = X(1,:) < obj.xgMin | X(2,:) < obj.ygMin | X(1,:) > obj.xgMax | X(2,:) > obj.ygMax;
      X(:,out_of_bounds) = [];
      X = sortrows(X',[3,1,2]);  % sorting is important for later handling
      
      %       % diagnostic
      %       plot(obj,'.b'); hold on;
      %       plot( [x x+dx]', [y y+dy]', 'x-b' )
      %       plot(X(:,1), X(:,2), 'xm');
      %       hold off
    end
    
    function plot(obj, varargin)
      %X = combinations( obj.gx, obj.gy );
      gridList = []; % array of [i,j] grid cell coordinates to plot
      for k=1:2:length(varargin)-1
        if strcmp( varargin{k}, 'gridList' )
          gridList = varargin{k+1};
        end
      end
      
      if isempty(gridList)
        X = obj.xg;  Y = obj.yg;
      else
        I = union( gridList(:,1), gridList(:,1) + 1 );
        J = union( gridList(:,2), gridList(:,2) + 1 );
        X = obj.xg(I);  Y = obj.yg(J);
      end
      
      % vertical lines
      plot( [X; X], [repmat(Y(1), 1, length(X)); repmat(Y(end), 1, length(X))], '-b', varargin{:} );
      hold on
      % horizontal lines
      plot( [repmat(X(1), 1, length(Y)); repmat(X(end), 1, length(Y))], [Y; Y], '-b', varargin{:} );
      hold off
    end
    
    function [imgdat, h] = imagesc2( obj, IJ, Z )
      imgdat.x = obj.xc;
      imgdat.y = obj.yc;
      imgdat.z = zeros( obj.Nx, obj.Ny );
      imgdat.z( sub2ind(size(imgdat.z), IJ(:,1), IJ(:,2)) ) = Z;
      imgdat.z = imgdat.z';
      h = imagesc( obj.xc, obj.yc, imgdat.z );
      set(h, 'AlphaData', imgdat.z~=0);
    end
    
    function [imgdat, h] = imagesc( obj, Z, alpha_value )
      imgdat.x = obj.xc;
      imgdat.y = obj.yc;
      if all(size(Z)==[obj.Nx obj.Ny])
        imgdat.z = Z';
      elseif all(size(Z)==[obj.Ny obj.Nx])
        imgdat.z = Z;
      else
        error('Z has invalid size (%d x %d, need %d x %d)', ...
          size(Z,1), size(Z,2), obj.Nx, obj.Ny);
      end
      h = imagesc( obj.xc, obj.yc, imgdat.z );
      if exist('alpha_value','var')
        set(h, 'AlphaData', imgdat.z ~= alpha_value);
      end
      axis xy; colorbar
    end
    
    function [i,j] = coord2sub( obj, x, y )
      % convert coordinates x, y to grid index i, j (starting with 1,1 in
      % the bottom left)
      % Return [NaN,NaN] if point [x,y] lies outside the grid.
      % Return half-integral numbers if point lies on an edge, e.g. [3.5,2] if
      % point lies on the edge between cell [3,2] and [4,2].
      i = (x-obj.xcMin+obj.ex/2)/obj.ex + 1;
      iEdge = i==round(i); i(iEdge) = i(iEdge)-0.5; i(~iEdge) = floor(i(~iEdge));
      j = (y-obj.ycMin+obj.ey/2)/obj.ey + 1;
      jEdge = j==round(j); j(jEdge) = j(jEdge)-0.5; j(~jEdge) = floor(j(~jEdge));
      out_of_bounds = x < obj.xgMin | y < obj.ygMin | x > obj.xgMax | y > obj.ygMax;
      i(out_of_bounds) = NaN;
      j(out_of_bounds) = NaN;
    end
    
    function ind = sub2ind( obj, I, J )
      ind = sub2ind( size(obj), I, J );
    end
    
    function ind = coord2ind( obj, x, y )
      %% convert coordinates to linear grid index
      [i,j] = obj.coord2sub( x, y );
      ind = obj.sub2ind( i, j );
    end
    
    function coord = ind2coord( obj, ind )
      [i,j]=ind2sub( size(obj), ind );
      coord = [obj.xc(i)' obj.yc(j)'];
    end
    
    function [S, O] = line_density( obj, lines )
      % lines: rows of [x y dx dy]
      % S: Nx-Ny-matrix of density inside grid
      % O: density outside grid
      
      % sort lines
      lines = normalize_lines( lines );
      % find intersections with grid
      X = obj.intersectLineSegments( lines );
      % translate to grid indices
      [xi, xj] = obj.coord2sub( X(:,1), X(:,2) );
      XI = [xi,xj];
      % S: cumulated line lengths (one entry per grid box)
      S = zeros(obj.Nx, obj.Ny);
      O = 0;
      for lineIndex = unique(X(:,3))'  % lines spanning multiple boxes
        line = lines(lineIndex,:);
        line_x0 = line(1:2);  % start point
        line_x1 = line_x0 + line(3:4);  % end point
        % intersection point indices with grid of this line
        I = X(:,3)==lineIndex;
        %yLine = Y(I,:);
        % number of intersections for this line
        n = nnz(I);
        % the relevant (for this line) intersection points
        % (should be sorted!)
        xi = XI(I,:); x = X(I,1:2);
        % for each line segment (i.e. for each box),
        % we calculate and add the line length
        assert( size(xi,1) == n );
        for ls = 1:(n+1)
          % find the box index and length of this line segment
          if ls==1  % one point at the beginning
            [i,j] = obj.coord2sub( line_x0(:,1), line_x0(:,2) );
            segment_length = norm( line_x0 - x(1,:) );
          elseif ls==n+1  % one point at the end
            [i,j] = obj.coord2sub( line_x1(:,1), line_x1(:,2) );
            segment_length = norm( line_x1 - x(end,:) );
          else  % two points somewhere in-between
            y = x(ls-1:ls,:);
            segment_length = norm( diff( x(ls-1:ls,1:2) ) );
            m = mean(y);
            [i,j] = obj.coord2sub(m(1),m(2));
          end
          if ~any(isnan([i,j]))
            if mod(i,1) > 0 || mod(j,1) > 0
              % we have a non-integer coordinate => try to skip
              assert( segment_length < 100*eps );
            else
              S(i,j) = S(i,j) + segment_length;
            end
          else
            O = O + segment_length;
          end
          % fprintf('Line #%d/Segment #%d at grid location (%d,%d).\n', lineIndex, ls, i, j);
        end
        %I = lines(:,:,lines(:,3)==lineIndex)
        
      end
      
      % these lines are contained completely in a box
      in_a_box = ~ismember(1:size(lines,1), unique(X(:,3)));
      for lineIndex = find(in_a_box)
        line = lines(lineIndex,:);
        line_x0 = line(1:2);  % start point
        [i,j] = obj.coord2sub( line_x0(:,1), line_x0(:,2) );
        segment_length = norm( line(3:4) );
        if ~any(isnan([i,j]))
          S(i,j) = S(i,j) + segment_length;
        else
          O = O + segment_length;
        end
        % fprintf('Line #%d contained in single grid location (%d,%d).\n', lineIndex, i, j);
      end
    end
    
    function [X,Y] = ndgrid( obj )
      [X,Y] = ndgrid(obj.xc, obj.yc);
    end
    
    function grid_width = grid_width( obj )
      assert( obj.ex == obj.ey, 'require quadratic grid');
      grid_width = obj.ex;
    end
    
    function X = coordinates( obj )
      %% return the grid coordinates as N-by-2 matrix.
      [x,y] = obj.ndgrid;
      X = [x(:) y(:)];
    end
    
    function bb = bounding_box( obj )
      %% bounding_box in format
      % [minX minY; maxX maxY]
      bb = [obj.xgMin obj.ygMin; obj.xgMax obj.ygMax];
    end
    
    function bbi = bounding_box_ind( obj )
      %% bounding_box_ind in format
      % [minXi minYi; maxXi maxYi]
      bbi = [1 1; obj.Nx obj.Ny];
    end
  end
  
end
