function fit = fit( obj, marked_network, varargin )
% AbstractModel#fit

ROWNAMEFORMAT = '%d: %d %.6f';

% parse options
argi = 1;
while argi <= numel(varargin)
  switch varargin{argi}
    case 'name'
      argi = argi + 1; name = varargin{argi};
    case 'wiggle_amount'
      argi = argi + 1; wiggle_amount = varargin{argi};
    case 'N_dummies'
      argi = argi + 1; N_dummies = varargin{argi};
    case 'N_simulations'
      argi = argi + 1; N_simulations = varargin{argi};
    case 'SaveImages'
      argi = argi + 1; SaveImages = varargin{argi};
    case 'SaveRealizationImages'
      argi = argi + 1; SaveRealizationImages = varargin{argi};
    case 'SaveState'
      argi = argi + 1; SaveState = varargin{argi};
    otherwise
      error('Unkown option %s',varargin{argi});
  end
  argi = argi + 1;
end

% set defaults
if isa( marked_network, 'MarkedNetwork' )
  if ~exist('N_dummies','var'); N_dummies = 3*marked_network.N; end
elseif isa( marked_network, 'function_handle' )
  marked_network_fun = marked_network;
  if ~exist('N_dummies','var'); error('please supply N_dummies'); end
else
  error('Invalid marked_network parameter.');
end
if ~exist('marked_network_fun','var'); marked_network_fun = []; end

if ~exist('name','var'); error('require a nice readable name for this fit'); end
if ~exist('wiggle_amount','var'); wiggle_amount = false; end
if ~exist('N_simulations','var'); N_simulations = 1; end
if ~exist('SaveImages','var'); SaveImages = false; end
if ~exist('SaveRealizationImages','var'); SaveRealizationImages = false; end
if ~exist('SaveState','var'); SaveState = false; end

configvars = {
  % config vars
  'N_cases'
  'wiggle_amount'
  'N_dummies'
  'N_simulations'
  'SaveImages'
  'SaveRealizationImages'
  };

% save variable images
if SaveImages && ismethod(obj,'save_variable_images'); obj.save_variable_images; end

% initialize network
if isa(obj,'Network.Model')
  net = obj.network;
else
  net = marked_network.net;
end

% find number of cases
if exist('marked_network_fun','var') && isa(marked_network_fun, 'function_handle')
  edgelist_marks = marked_network_fun();
  N_cases = size(edgelist_marks,1);
else
  N_cases = marked_network.nummarks;
end

% calculate independent variables at mark locations (data & dummies)
varnames = {obj.variables.name};
variables = obj.variables;

%% load/compute model realizations
if exist( filename(name,'variates.mat'), 'file' ) && ...
    exist( filename(name,'realizations.mat'), 'file' )
  fprintf('Loading variates & realizations from file...\n');
  load( filename(name,'variates.mat') );
  load( filename(name,'realizations.mat') );
end

if exist('variates','var') && exist('realizations','var')
  sim0 = min(size(variates,1), size(realizations,1)) + 1;
else
  sim0 = 1;
  realizations = [];
  variates = struct([]);
end

if N_simulations >= sim0
  nNewSims = N_simulations-sim0+1;
  fprintf('Creating %d model realizations: ', nNewSims);

  realizations = [realizations; cell(nNewSims,1)];

  empty_variate_struct = struct('W',[], 'Xedge',[], 'YW',[]);
  variates = [variates; ...
    repelem( empty_variate_struct, nNewSims, 1 )];

  assert( numel(realizations) == N_simulations );
  assert( size(variates,1) == N_simulations );

  parfor sim = sim0:N_simulations

    if isa(marked_network_fun, 'function_handle')
      % we have a function to generate marked networks on-the-fly
      edgelist_marks = marked_network_fun();
    else
      % wiggle cases (recommended for numerical stability)
      if wiggle_amount
        edgelist_marks = net.rmove( marked_network.edge_network, wiggle_amount );
      else
        edgelist_marks = marked_network.edge_network;
      end
    end

    % generate dummy points
    edgelist_dummies = net.rpois(N_dummies);

    % fix coordinate ordering
    Y = [ones(N_cases,1); zeros(N_dummies,1)];

    % build combined marked network of both data & dummies
    mn = MarkedNetwork( net, [edgelist_marks; edgelist_dummies] );

    % assert simpleness
    assert( mn.is_simple, 'Network is not simple. Cannot continue. Use wiggle_amount > 0 to fix.' );

    % build independent and dependent variable vector/matrix
    Y = Y(mn.edge_coordinates_ordering,:);

    % calculate weights
    variates(sim).W = mn.weights;
    if isa(variates(sim).W,'NetworkTessellation')
      variates(sim).W = variates(sim).W.weights;
    end

    % check weights
    if min(variates(sim).W)>0 && max(variates(sim).W)~= Inf
    else
      warning('Retrying because of invalid tesselation...');
      continue  % retry
    end

    % optionally, save model realization images
    if SaveRealizationImages && isa(obj,'Grid.Model')
      Coords = net.xy_from_edge( edgelist_marks(:,1), edgelist_marks(:,2) );
      [~,~,~,CountsMatrix] = obj.G.boxcount( Coords(:,1), Coords(:,2) );
      obj.G.save_image( CountsMatrix, filename(name,'model_realization'), 'partial', sim );
    end

    % save model realization to var
    realizations{sim} = edgelist_marks;

    % find edgelist format of marks
    variates(sim).Xedge = EdgeCoordinates( mn.net, mn.edge_network );

    % build dependent variable
    variates(sim).YW = Y ./ variates(sim).W;

    % print progress
    fprintf('.\n');
  end
  % save intermediate results
  %variates = cell2mat( variates );
  if SaveState
    save(filename(name,'realizations.mat'), 'realizations', '-v7.3');
    save(filename(name,'variates.mat'), 'variates', '-v7.3');
  end
  fprintf('\n');
end

grpSize = numel(variates(1).YW);

if true

  %% load/compute covariates
  if exist( filename(name,'flat-covariates.mat'), 'file' )
    fprintf('Loading covariates from flat file...\n');
    load( filename(name,'flat-covariates.mat') );
  end

  Xedges = EdgeCoordinates( net, ...
    cell2mat( arrayfun(@(v) v.Xedge.edge_network, ...
    variates, 'UniformOutput', 0) ) );

  if ~exist('TX','var') || ~istable(TX) || table_overvations_missing(TX)
    TX = table;
  end

  if table_variables_missing(TX)
    fprintf('Building (missing) covariates: ');
    [~,missingIdx] = table_variables_missing(TX);
    colindex = size(TX,2) + 1;
    for model_index = missingIdx
      col = obj.edge_value( variables(model_index), Xedges );
      TX.(colindex) = col;
      TX.Properties.VariableNames{colindex} = varnames{model_index};
      fprintf('%s complete.\n', varnames{model_index});
      colindex = colindex + 1;
    end

    % save intermediate results
    if SaveState
      save(filename(name,'flat-covariates.mat'), 'TX', '-v7.3');
    end
  end

  fprintf('\nValidate covariates...\n');
  for model_index = 1:size(TX,2)
    % validate independent variables
    col = TX.(model_index);
    if ~iscategorical(col)
      if all( col == 0 ); error('all-zero in %s', varnames{model_index}); end
      if any( isinf(col) ); varerror('Inf', varnames{model_index}, find(isinf(col),1)); end
      if any( isnan(col) ); varerror('NaN', varnames{model_index}, find(isnan(col),1)); end
    end
  end
  % assert_correct_rownames( TX, Xedges );

end

fprintf('Reorganize covariates: ');
for sim = 1:N_simulations
  variates(sim).TX = TX( (sim-1)*grpSize + (1:grpSize), : );
  fprintf('.');
end
clear TX
fprintf('\n');

fprintf('Perform fitting: ');
results = cell(N_simulations, 1);
for sim = 1:N_simulations

  % fit model
  T = variates(sim).TX(:,varnames);
  T.YW = variates(sim).YW;
  model = fitglm( T, ...
    'PredictorVars', varnames, ...
    'CategoricalVars', varnames([variables.categorical]), ...
    'Distribution', 'poisson', ...
    'DispersionFlag', true, 'Weights', variates(sim).W, ...
    'Intercept', false, 'Offset', 0, 'Link', 'log' );

  % save results
  theta = model.Coefficients.Estimate';
  results{sim} = struct('theta',theta, 'meta',[model.Dispersion],...
    'meta_names', {'dispersion'}, 'fitglm_output', model);

  % log progress
  fprintf('%3d: dispersion = %g', sim, model.Dispersion);
  for j = 1:numel(model.CoefficientNames)
    fprintf(', %s = %.2e', model.CoefficientNames{j}, theta(j));
  end
  fprintf('\n');

end
% save intermediate results
results = cell2mat(results);
if SaveState
  save(filename(name,'results.mat'), 'results', '-v7.3');
end


state = saveState;

% return ModelFit object
switch class(obj)
  case 'Grid.Model'
    fit = Grid.ModelFit( name, obj, results, state, realizations, variates );
  case 'Network.Model'
    % FIXME: this will fail / marked_network can be a function_handle
    fit = Network.ModelFit( marked_network, name, obj, results, state, realizations, variates );
  otherwise
    error('Unknown Model class: %s',class(obj));
end


  function state = saveState()
    % http://stackoverflow.com/a/3470731
    values = cellfun(@(n) evalin('caller',n),configvars,'UniformOutput',false);
    state = cell2struct(values,configvars,1);
    %if SaveState
    %  save(StateFile, configvars{:}, '-v7.3');
    %end
  end

  function varerror( msg, varname, erridx )
    xy = net.xy_from_edge( Xedges.edge_network(erridx,1), ...
      Xedges.edge_network(erridx,2) );
    warning('%s in %s [#%d, %f x %f]', msg, ...
      varname, erridx, xy(1), xy(2) );
  end

  function [tf, missingIdx] = table_variables_missing( T )
    % check completeness of variables
    varPresent = ismember( varnames, T.Properties.VariableNames);
    tf = ~all( varPresent );
    if nargout > 1
      missingIdx = find(~varPresent);
    end
  end

  function tf = table_overvations_missing( T )
    % check completeness of observations
    tf = size(T,1) < N_simulations * grpSize;
  end

  function cellstr = XedgeToString
    cellstr = arrayfun( @(i) sprintf(ROWNAMEFORMAT, Xedges.edge_network(i,1), Xedges.edge_network(i,2)), 1:Xedges.numpoints, 'UniformOutput', false );
  end

  function assert_correct_rownames( T, Xedges )
    % check for correctness of observations
    assert( numel(T.Properties.RowNames) == Xedges.numpoints );
    assert( all(strcmp(XedgeToString', T.Properties.RowNames)) );
  end

end
