function filename = filename( varargin )

filename = strrep( strjoin(varargin,'-'), ' ', '_' );
