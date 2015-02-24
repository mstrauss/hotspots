function intensity = intensity( obj, save )

assert( isa(obj.model,'Grid.Model') );

G = obj.model.G;
Gvar = obj.model.Gvariables;
theta = obj.theta;

% the names of the variables where we have coefficients
varnames = obj.varnames;
assert( numel(varnames) == numel(theta) );

% the names of the variables where we have covariate values
covarnames = obj.model.Gvarnames;
assert( numel(covarnames) == numel(Gvar) );

[~,covarindex] = ismember( varnames, covarnames );
assert( all( covarindex > 0 ), 'not all fit variables present in model');
    

if ~exist('save','var'); save = false; end

fprintf('Calculating intensity for a %s model using variables/coefficients:\n  %s\n', ...
  obj.model.type, strjoin( cellfun( @(v,c) sprintf('%s = %g',v,c), varnames, num2cell(theta), 'UniformOutput', 0 ), ', ') );


switch obj.model.type
  case 'loglinear'
    % init.
    intensity = zeros(size(G));

    % the network length (per pixel)
    L = obj.model.G.line_density(obj.model.network.edges);  
    
    % loop over fit variable (for each categorical var in the model,
    % there is one fit var. for each level)
    for fit_index = 1:numel(varnames)
      gi = covarindex(fit_index);
      varvalue = Gvar(gi).data;
      if Gvar(gi).variable.categorical
        if isa(Gvar(gi).variable.domain,'Domain.NetworkEdge')
          % we use the ratio of the network lengths per pixel of the fitted
          % categorical variable vs. the total network lengt per pixel
          varvalue = varvalue ./ L;
        else
          error('unsupported type of categorical variable')
        end
      end
      intensity = intensity + varvalue * theta(fit_index);
    end
    
    if save
      %% plot model intensity
      obj.model.G.save_image( intensity, filename(obj.name, 'log-intensity without network'), 'grid_with', obj.model.G.ex );
      obj.model.G.save_image( L, filename(obj.name, 'network length'), 'grid_with', obj.model.G.ex );
    end

    if any(isinf(intensity(:))); warning('Inf in resulting intensity. Probably not what you want.'); end
    if any(isnan(intensity(:))); warning('NaN in resulting intensity. Probably not what you want.'); end
    % clean Inf and NaN
    intensity(isnan(intensity)) = -Inf;
    intensity(intensity==0) = -Inf;
    
    intensity = exp(intensity);
    
    if save
      obj.model.G.save_image( intensity, filename(obj.name, 'intensity without network'), 'grid_with', obj.model.G.ex );
    end

    intensity = intensity .* L;  % we live on the network domain

  case 'linear2'
    % linear model
    warning('THIS IS OBSOLETE');
    L = obj.model.G.line_density(obj.marked_network.net.edges);
    intensity = repmat( theta(1), G.Nx, G.Ny ) + ...
      GX{2} .* theta(2) + ...
      GX{3} .* theta(3);
    intensity = intensity .* GX{4} * theta(4);
    %   for i = numel(varnames)
    %     model_intensity = model_intensity + ...
    %       GX(:,:,i) .* T.C(varnames(i));
    %   end
  otherwise
    error( 'Unknown model type: %s', obj.model.type );
end

if save
  %% plot model intensity
  obj.model.G.save_image( intensity, filename(obj.name, 'intensity') );
  obj.model.G.save_image( log(intensity), filename(obj.name, 'log-intensity') );
end
