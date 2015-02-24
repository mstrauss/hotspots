function model_intensity = intensity( obj )
%% Network.ModelFit#intensity calculates the intensity function on the network.
% First, we calculate a network tesselation and measure the covariates
% there;  On the cells of this network tesselation, we compute the
% intensities from the covariates and the cell length.  The quality of the
% tesselation (and the assumed "constantness" of the covariates within a
% single cell) shall be computed a-posteriori from a TBD smoothness
% criterion.

theta = obj.theta;
net = obj.model.network;
mean_resolution = 10;  % (km) mean linear resolution
nCells = round(net.length / mean_resolution);

switch obj.model.type
  case 'loglinear'
    %% log-linear model
    
    % generate tesselation
    [~,mn] = net.rpois(nCells);  % FIXME: need better tesselation algo.
    
    % assert simpleness
    assert( mn.is_simple );
    
    % calculate weights
    [W,U] = mn.weights;
    assert( abs( U + sum(W) - net.length )/net.length < 1e-10 );
    
    %% find covariates at cell center (=marks) positions
    nVars = numel(obj.varnames);
    X = zeros(nCells, nVars );
    Xedges = EdgeCoordinates( net, mn.edge_network );
    fit_index = 1;
    for model_index = 1:numel(obj.model.variables)  % loop over model variables
      mv = obj.model.variables(model_index);  % model variable
      val = obj.model.variables(model_index).edge_value( Xedges );  % variable value
      if mv.categorical
        valcat = categories(val);
        for cati = 1:numel(valcat)
          cat = valcat(cati);
          if ~ismember(sprintf('%s_%d', mv.name, cati), obj.varnames); continue; end
          X(:,fit_index) = val==cat;
          fit_index = fit_index + 1;
        end
      else
        % use value 1:1
        X(:,fit_index) = val;
        fit_index = fit_index + 1;
      end
      fprintf('%s complete.\n', obj.model.varnames{model_index});
    end
    
    % make table
    TX = array2table(X,'VariableNames',obj.varnames);
    
    % add cell weights & coordinates
    TX.W = W;
    TX.XY = mn.coords;
    TX.ED = mn.edge_network;
    
    TX.lambda = zeros(nCells,1);
    for vi = 1:nVars
      TX.lambda = TX.lambda + theta(vi) * TX.(vi);
    end
    TX.prob = exp(TX.lambda) .* TX.W;
    TX.prob = TX.prob / sum(TX.prob);
    
    model_intensity = ModelIntensity( Xedges, TX.prob );
    
  case 'linear2'
    % linear model
    error('THIS IS OBSOLETE');
  otherwise
    error( 'Unknown model type: %s', obj.model.type );
end
